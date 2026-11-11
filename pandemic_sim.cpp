#include <thread>
#include <future>
#include <queue>
#include <mutex>
#include <condition_variable>
#include <functional>
#include <memory>
#include <type_traits>
#include <stdexcept>
#include <utility>
#include <iostream>
#include <vector>
#include <string>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <numeric>
#include <algorithm>
#include <random>

// State definitions: S, I, R, V, D
enum State : uint8_t {
    STATE_S = 0,
    STATE_I = 1,
    STATE_R = 2,
    STATE_V = 3,
    STATE_D = 4
};

struct Cell {
    uint8_t state;
    uint8_t c; // days infected counter
};

struct Counts {
    long long S = 0;
    long long I = 0;
    long long R = 0;
    long long V = 0;
    long long D = 0;
};

struct Region {
    int rowStart;
    int rowEnd;
    int colStart;
    int colEnd;
    int rowStride = 1;
};

// Simulation Parameters
struct Params {
    int N = 1000;               // Grid size N x N
    int max_days = 200;         // Maximum number of days
    double beta = 0.08;         // Transmission probability
    int d = 7;                  // Days infectious
    double mu = 0.01;           // Mortality
    int I0 = 10;                // Initially infected
    double v = 0.0;             // Vaccinated fraction
    double m = 1.0;             // Mask factor (0 < m <= 1)
    double lambda = 0.3;        // Lockdown factor
    double theta_on = 0.02;     // Lockdown threshold ON
    double theta_off = 0.005;   // Lockdown threshold OFF
    uint64_t seed = 42;         // Random seed
    int num_threads = 1;        // Number of threads T (Part A sequential = 1)
    int num_chunks = 1;         // Number of chunks C
    int snapshot_interval = 0;  // Snapshot interval K (0 = off)
    int snapshot_limit = 4;     // Limit Q pending writes
    int io_latency_ms = 0;      // --io-latency in ms
    bool fixed_days = false;    // Whether to force running fixed number of days
    std::string output_file = "output.csv"; // CSV output path
    bool use_thread_pool = false;
    std::string partition = "strips";
};


class ThreadPool {
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    bool stopping = false;

public:
    explicit ThreadPool(size_t threadCount) {
        if (threadCount == 0) {
            throw std::invalid_argument("Thread count must be positive");
        }

        for (size_t i = 0; i < threadCount; ++i) {

            workers.emplace_back([this]() {

                while (true) {
                    std::function<void()> task;

                    {
                        std::unique_lock<std::mutex> lock(queueMutex);

                        condition.wait(lock, [this]() {
                            return stopping || !tasks.empty();
                        });

                        if (stopping && tasks.empty()) {
                            return;
                        }

                        task = std::move(tasks.front());
                        tasks.pop();
                    }

                    task();
                }
            });
        }
    }

    template <typename F>
    auto submit(F&& function)
        -> std::future<std::invoke_result_t<F>> {

        using ReturnType = std::invoke_result_t<F>;

        auto task = std::make_shared<
            std::packaged_task<ReturnType()>
        >(std::forward<F>(function));

        std::future<ReturnType> result = task->get_future();

        {
            std::lock_guard<std::mutex> lock(queueMutex);

            if (stopping) {
                throw std::runtime_error("Thread pool is stopping");
            }

            tasks.emplace([task]() {
                (*task)();
            });
        }

        condition.notify_one();

        return result;
    }

    ~ThreadPool() {

        {
            std::lock_guard<std::mutex> lock(queueMutex);
            stopping = true;
        }

        condition.notify_all();

        for (auto& worker : workers) {
            if (worker.joinable()) {
                worker.join();
            }
        }
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;
};


// High-quality hash-based PRNG returning uniform double in [0, 1)
// Function of (seed, day, cell, stream) as specified in Section 3.1
inline uint64_t hash_coords(uint64_t seed, uint64_t day, uint64_t cell, uint64_t stream) {
    uint64_t h = seed ^ (day * 0x9e3779b97f4a7c15ULL);
    h ^= (cell * 0xbf58476d1ce4e5b9ULL) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    h ^= (stream * 0x94d049bb133111ebULL) + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    // SplitMix64 finalizer
    h ^= h >> 30;
    h *= 0xbf58476d1ce4e5b9ULL;
    h ^= h >> 27;
    h *= 0x94d049bb133111ebULL;
    h ^= h >> 31;
    return h;
}

inline double rand01(uint64_t seed, uint64_t day, uint64_t cell, uint64_t stream) {
    uint64_t h = hash_coords(seed, day, cell, stream);
    return (h >> 11) * (1.0 / 9007199254740992.0); // 1.0 / 2^53
}

// 64-bit FNV-1a hash over the entire grid
uint64_t computeGridChecksum(const std::vector<Cell>& grid) {
    uint64_t hash = 14695981039346656037ULL;
    for (const auto& cell : grid) {
        uint16_t val = (static_cast<uint16_t>(cell.state) << 8) | cell.c;
        hash ^= (val & 0xFF);
        hash *= 1099511628211ULL;
        hash ^= ((val >> 8) & 0xFF);
        hash *= 1099511628211ULL;
    }
    return hash;
}



Counts updateRegion(
    const std::vector<Cell>& cur,
    std::vector<Cell>& next,
    const Region& region,
    int N,
    int day,
    const double P[9],
    const Params& p
) {
    Counts localCounts;

    const int di[8] = {-1,-1,-1,0,0,1,1,1};
    const int dj[8] = {-1,0,1,-1,1,-1,0,1};

    for (int i = region.rowStart;
     i < region.rowEnd;
     i += region.rowStride){
        for (int j = region.colStart; j < region.colEnd; ++j) {

            long long id = 1LL * i * N + j;
            const Cell& cell = cur[id];

            if (cell.state == STATE_S) {

                int k = 0;

                for (int dir = 0; dir < 8; ++dir) {
                    int ni = i + di[dir];
                    int nj = j + dj[dir];

                    if (ni >= 0 && ni < N &&
                        nj >= 0 && nj < N &&
                        cur[1LL * ni * N + nj].state == STATE_I) {
                        k++;
                    }
                }

                if (k > 0 &&
                    rand01(p.seed, day, id, 0) < P[k]) {
                    next[id] = {STATE_I, 0};
                    localCounts.I++;
                } else {
                    next[id] = {STATE_S, 0};
                    localCounts.S++;
                }

            } else if (cell.state == STATE_I) {

                if (cell.c + 1 == p.d) {

                    if (rand01(p.seed, day, id, 1) < p.mu) {
                        next[id] = {STATE_D, 0};
                        localCounts.D++;
                    } else {
                        next[id] = {STATE_R, 0};
                        localCounts.R++;
                    }

                } else {
                    next[id] = {
                        STATE_I,
                        static_cast<uint8_t>(cell.c + 1)
                    };
                    localCounts.I++;
                }

            } else {

                next[id] = cell;

                if (cell.state == STATE_R) localCounts.R++;
                else if (cell.state == STATE_V) localCounts.V++;
                else if (cell.state == STATE_D) localCounts.D++;
            }
        }
    }

    return localCounts;
}


// Keep this wrapper for existing B1 and B2 code
Counts updateChunk(
    const std::vector<Cell>& cur,
    std::vector<Cell>& next,
    int startRow,
    int endRow,
    int N,
    int day,
    const double P[9],
    const Params& p
) {
    Region region{startRow, endRow, 0, N};

    return updateRegion(
        cur, next, region, N, day, P, p
    );
}

void printHelp(const char* prog) {
    std::cout << "Usage: " << prog << " [options]\n"
              << "Options:\n"
              << "  -N, --grid-size <int>         Grid dimension N (default: 1000)\n"
              << "  -D, --days <int>              Number of days to simulate (default: 200)\n"
              << "  -s, --seed <uint64>           Random seed (default: 42)\n"
              << "  -b, --beta <float>            Transmission probability beta (default: 0.08)\n"
              << "  -d, --days-infectious <int>   Days infectious d (default: 7)\n"
              << "  -u, --mu, --mortality <float> Mortality probability mu (default: 0.01)\n"
              << "  --I0, --initially-infected <int> Initially infected count (default: 10)\n"
              << "  -v, --vaccinated <float>      Vaccinated fraction v (default: 0.0)\n"
              << "  -m, --mask-factor <float>     Mask factor m (default: 1.0)\n"
              << "  -l, --lambda <float>          Lockdown factor lambda (default: 0.3)\n"
              << "  --theta-on <float>            Lockdown threshold ON (default: 0.02)\n"
              << "  --theta-off <float>           Lockdown threshold OFF (default: 0.005)\n"
              << "  -T, --threads <int>           Number of threads (default: 1)\n"
              << "  -C, --chunks <int>            Number of chunks (default: 1)\n"
              << "  -K, --snapshot-interval <int> Snapshot interval in days (default: 0)\n"
              << "  -Q, --snapshot-limit <int>    Max pending snapshot writes (default: 4)\n"
              << "  --io-latency <int>            Simulated I/O latency in ms (default: 0)\n"
              << "  --fixed-days                  Simulate exact number of days even if I reaches 0\n"
              << "  -o, --output <file>           Output CSV filename (default: output.csv)\n"
              << "  -h, --help                    Show this help message\n"
              << "  --pool                       Use reusable thread pool\n"
              << "  --partition <type>    strips, blocks or interleaved\n";
}

Params parseArgs(int argc, char* argv[]) {
    Params p;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if ((arg == "-N" || arg == "--grid-size") && i + 1 < argc) {
            p.N = std::stoi(argv[++i]);
        } else if ((arg == "-D" || arg == "--days") && i + 1 < argc) {
            p.max_days = std::stoi(argv[++i]);
        } else if ((arg == "-s" || arg == "--seed") && i + 1 < argc) {
            p.seed = std::stoull(argv[++i]);
        } else if ((arg == "-b" || arg == "--beta") && i + 1 < argc) {
            p.beta = std::stod(argv[++i]);
        } else if ((arg == "-d" || arg == "--days-infectious") && i + 1 < argc) {
            p.d = std::stoi(argv[++i]);
        } else if ((arg == "-u" || arg == "--mu" || arg == "--mortality") && i + 1 < argc) {
            p.mu = std::stod(argv[++i]);
        } else if ((arg == "--I0" || arg == "--initially-infected") && i + 1 < argc) {
            p.I0 = std::stoi(argv[++i]);
        } else if ((arg == "-v" || arg == "--vaccinated") && i + 1 < argc) {
            p.v = std::stod(argv[++i]);
        } else if ((arg == "-m" || arg == "--mask-factor") && i + 1 < argc) {
            p.m = std::stod(argv[++i]);
        } else if ((arg == "-l" || arg == "--lambda") && i + 1 < argc) {
            p.lambda = std::stod(argv[++i]);
        } else if (arg == "--theta-on" && i + 1 < argc) {
            p.theta_on = std::stod(argv[++i]);
        } else if (arg == "--theta-off" && i + 1 < argc) {
            p.theta_off = std::stod(argv[++i]);
        } else if ((arg == "-T" || arg == "--threads") && i + 1 < argc) {
            p.num_threads = std::stoi(argv[++i]);
        } else if ((arg == "-C" || arg == "--chunks") && i + 1 < argc) {
            p.num_chunks = std::stoi(argv[++i]);
        } else if ((arg == "-K" || arg == "--snapshot-interval") && i + 1 < argc) {
            p.snapshot_interval = std::stoi(argv[++i]);
        } else if ((arg == "-Q" || arg == "--snapshot-limit") && i + 1 < argc) {
            p.snapshot_limit = std::stoi(argv[++i]);
        } else if (arg == "--io-latency" && i + 1 < argc) {
            p.io_latency_ms = std::stoi(argv[++i]);
        } else if (arg == "--fixed-days") {
            p.fixed_days = true;
        } else if ((arg == "-o" || arg == "--output") && i + 1 < argc) {
            p.output_file = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            printHelp(argv[0]);
            std::exit(0);
        } else if (arg == "--pool") {
            p.use_thread_pool = true;
        } else if (arg == "--partition" && i + 1 < argc) {
            p.partition = argv[++i];
        } else {
            std::cerr << "Warning: Unknown option " << arg << "\n";
        }
    }
    return p;
}

int main(int argc, char* argv[]) {
    Params p = parseArgs(argc, argv);

    const int N = p.N;
    const long long total_cells = static_cast<long long>(N) * N;

    // Grid buffers: current and next
    std::vector<Cell> cur(total_cells, {STATE_S, 0});
    std::vector<Cell> next(total_cells, {STATE_S, 0});

    // Initial conditions setup (done sequentially)
    std::mt19937_64 rng(p.seed);

    // 1. Vaccinate fraction v
    long long num_vaccinated = static_cast<long long>(std::round(p.v * total_cells));
    if (num_vaccinated > total_cells) num_vaccinated = total_cells;

    if (num_vaccinated > 0) {
        std::vector<long long> all_indices(total_cells);
        std::iota(all_indices.begin(), all_indices.end(), 0LL);
        // Partial Fisher-Yates shuffle
        for (long long i = 0; i < num_vaccinated; ++i) {
            std::uniform_int_distribution<long long> dist(i, total_cells - 1);
            long long j = dist(rng);
            std::swap(all_indices[i], all_indices[j]);
            cur[all_indices[i]].state = STATE_V;
        }
    }

    // 2. Infect I0 susceptible individuals chosen at random
    long long susceptible_count = total_cells - num_vaccinated;
    long long to_infect = std::min(static_cast<long long>(p.I0), susceptible_count);
    long long infected_so_far = 0;

    while (infected_so_far < to_infect) {
        std::uniform_int_distribution<long long> dist(0, total_cells - 1);
        long long idx = dist(rng);
        if (cur[idx].state == STATE_S) {
            cur[idx].state = STATE_I;
            cur[idx].c = 0;
            infected_so_far++;
        }
    }

    // Open output CSV
    std::ofstream csv(p.output_file);
    if (!csv.is_open()) {
        std::cerr << "Error: Could not open output file " << p.output_file << "\n";
        return 1;
    }
    csv << "day,S,I,R,V,D,beta_eff,checksum\n";

    // Initial counts for day 0
    Counts counts;
    for (long long idx = 0; idx < total_cells; ++idx) {
        switch (cur[idx].state) {
            case STATE_S: counts.S++; break;
            case STATE_I: counts.I++; break;
            case STATE_R: counts.R++; break;
            case STATE_V: counts.V++; break;
            case STATE_D: counts.D++; break;
        }
    }

    bool lockdown = false;
    double beta_eff = p.beta * p.m;
    uint64_t checksum = computeGridChecksum(cur);

    // Write Day 0 initial state
    csv << "0,"
        << counts.S << ","
        << counts.I << ","
        << counts.R << ","
        << counts.V << ","
        << counts.D << ","
        << std::fixed << std::setprecision(6) << beta_eff << ","
        << "0x" << std::hex << std::setw(16) << std::setfill('0') << checksum << std::dec << "\n";

    // Precompute P[k] table for current beta_eff
    // Moore neighborhood offsets

    std::unique_ptr<ThreadPool> pool;

    if (p.use_thread_pool) {
        pool = std::make_unique<ThreadPool>(p.num_threads);
    }

    // Main simulation loop
    for (int day = 0; day < p.max_days; ++day) {
        // Check stopping condition: stop if no infected and not fixed_days
        if (!p.fixed_days && counts.I == 0) {
            break;
        }

        // 1. Determine beta_eff for the upcoming day
        // Decisions made from day t's totals apply to day t + 1
        double infected_fraction = static_cast<double>(counts.I) / static_cast<double>(total_cells);
        if (!lockdown) {
            if (infected_fraction > p.theta_on) {
                lockdown = true;
            }
        } else {
            if (infected_fraction < p.theta_off) {
                lockdown = false;
            }
        }

        beta_eff = (p.beta * p.m) * (lockdown ? p.lambda : 1.0);

        // Precompute P[k] = 1 - (1 - beta_eff)^k for k = 0..8
        double P[9];
        for (int k = 0; k <= 8; ++k) {
            P[k] = 1.0 - std::pow(1.0 - beta_eff, k);
        }

        // Reset next counts
        Counts next_counts;


int C = std::max(1, p.num_chunks);
C = std::min(C, N);

std::vector<Region> regions;
regions.reserve(C);

if (p.partition == "strips") {

    // B3.1: Row strips
    int rowsPerChunk = N / C;
    int extraRows = N % C;
    int startRow = 0;

    for (int chunk = 0; chunk < C; ++chunk) {
        int rows = rowsPerChunk +
                   (chunk < extraRows ? 1 : 0);

        int endRow = startRow + rows;

        regions.push_back({startRow, endRow, 0, N});

        startRow = endRow;
    }

} else if (p.partition == "blocks") {

    // B3.2: Approximately square blocks
    int tileRows = static_cast<int>(std::ceil(std::sqrt(C)));

    int baseBlocks = C / tileRows;
    int extraBlocks = C % tileRows;

    int rowStart = 0;

    for (int r = 0; r < tileRows; ++r) {

        int rows = N / tileRows +
                   (r < N % tileRows ? 1 : 0);

        int rowEnd = rowStart + rows;

        int blocksInRow = baseBlocks +
                          (r < extraBlocks ? 1 : 0);

        int colStart = 0;

        for (int c = 0; c < blocksInRow; ++c) {

            int cols = N / blocksInRow +
                       (c < N % blocksInRow ? 1 : 0);

            int colEnd = colStart + cols;

            regions.push_back({
                rowStart, rowEnd,
                colStart, colEnd
            });

            colStart = colEnd;
        }

        rowStart = rowEnd;
    }

} else if (p.partition == "interleaved") {

    // B3.3: Interleaved rows
    int T = p.num_threads;

    if (T <= 0 || C < T) {
        throw std::invalid_argument(
            "Interleaved partition requires C >= T > 0"
        );
    }

    for (int chunk = 0; chunk < C; ++chunk) {

        int lane = chunk % T;
        int subtask = chunk / T;

        int tasksForLane = 1 + (C - 1 - lane) / T;

        int rowStart = lane + subtask * T;
        int stride = T * tasksForLane;

        regions.push_back({
            rowStart,
            N,
            0,
            N,
            stride
        });
    }

} else {
    throw std::invalid_argument(
        "Unknown partition. Use strips, blocks or interleaved."
    );
}



// Execute generated regions
next_counts = Counts{};

if (pool) {

    // B2: Thread pool
    std::vector<std::future<Counts>> futures;
    futures.reserve(regions.size());

    for (size_t chunk = 0; chunk < regions.size(); ++chunk) {

        futures.emplace_back(
            pool->submit(
                [&, chunk]() {
                    return updateRegion(
                        cur, next, regions[chunk],
                        N, day, P, p
                    );
                }
            )
        );
    }

    for (auto& future : futures) {
        Counts local = future.get();

        next_counts.S += local.S;
        next_counts.I += local.I;
        next_counts.R += local.R;
        next_counts.V += local.V;
        next_counts.D += local.D;
    }

} else {

    // B1: Thread per day
    std::vector<std::thread> workers;
    std::vector<Counts> threadCounts(regions.size());

    for (size_t chunk = 0; chunk < regions.size(); ++chunk) {

        workers.emplace_back(
            [&, chunk]() {
                threadCounts[chunk] = updateRegion(
                    cur, next, regions[chunk],
                    N, day, P, p
                );
            }
        );
    }

    for (auto& worker : workers) {
        worker.join();
    }

    for (const auto& local : threadCounts) {
        next_counts.S += local.S;
        next_counts.I += local.I;
        next_counts.R += local.R;
        next_counts.V += local.V;
        next_counts.D += local.D;
    }
}
        

        // Swap grids
        std::swap(cur, next);
        counts = next_counts;

        // Compute checksum of the newly updated grid
        checksum = computeGridChecksum(cur);

        // Write row for day + 1
        csv << (day + 1) << ","
            << counts.S << ","
            << counts.I << ","
            << counts.R << ","
            << counts.V << ","
            << counts.D << ","
            << std::fixed << std::setprecision(6) << beta_eff << ","
            << "0x" << std::hex << std::setw(16) << std::setfill('0') << checksum << std::dec << "\n";
    }

    csv.close();
    std::cout << "Simulation completed. Results written to " << p.output_file << "\n";
    return 0;
}

