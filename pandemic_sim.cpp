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
#include <chrono>
#include <deque>
#include <filesystem>

// State definitions: S, I, R, V, D
enum State : uint8_t
{
    STATE_S = 0,
    STATE_I = 1,
    STATE_R = 2,
    STATE_V = 3,
    STATE_D = 4
};

struct Cell
{
    uint8_t state;
    uint8_t c; // days infected counter
};

struct Counts
{
    long long S = 0;
    long long I = 0;
    long long R = 0;
    long long V = 0;
    long long D = 0;
};

struct Region
{
    int rowStart;
    int rowEnd;
    int colStart;
    int colEnd;
    int rowStride = 1;
};

// Simulation Parameters
struct Params
{
    int N = 1000;                           // Grid size N x N
    int max_days = 200;                     // Maximum number of days
    double beta = 0.08;                     // Transmission probability
    int d = 7;                              // Days infectious
    double mu = 0.01;                       // Mortality
    int I0 = 10;                            // Initially infected
    double v = 0.0;                         // Vaccinated fraction
    double m = 1.0;                         // Mask factor (0 < m <= 1)
    double lambda = 0.3;                    // Lockdown factor
    double theta_on = 0.02;                 // Lockdown threshold ON
    double theta_off = 0.005;               // Lockdown threshold OFF
    uint64_t seed = 42;                     // Random seed
    int num_threads = 1;                    // Number of threads T (Part A sequential = 1)
    int num_chunks = 1;                     // Number of chunks C
    int snapshot_interval = 0;              // Snapshot interval K (0 = off)
    int snapshot_limit = 4;                 // Limit Q pending writes
    int io_latency_ms = 0;                  // --io-latency in ms
    bool fixed_days = false;                // Whether to force running fixed number of days
    std::string output_file = "output.csv"; // CSV output path
    bool use_thread_pool = false;
    std::string partition = "strips";
    bool run_sweep = false;

    std::vector<double> sweep_betas = {0.05, 0.08, 0.10};
    std::vector<double> sweep_vaccinations = {0.0, 0.2, 0.4};
    std::vector<uint64_t> sweep_seeds = {42, 43, 44};

    int sweep_workers = 4;
    int time_limit_seconds = 0;

    std::string sweep_output =
        "results/part_c/parameter_sweep.csv";
    bool async_snapshots = false;
    std::string snapshot_prefix = "results/part_c/snapshot";
};

struct SweepSpec
{
    int run_id;
    double beta;
    double vaccination;
    uint64_t seed;
};

struct SweepSummary
{
    int run_id;
    double beta;
    double vaccination;
    uint64_t seed;

    double peak_infected_fraction;
    int day_of_peak;

    long long total_infected;
    long long total_deaths;

    int days;
};

class ThreadPool
{
private:
    std::vector<std::thread> workers;
    std::queue<std::function<void()>> tasks;

    std::mutex queueMutex;
    std::condition_variable condition;
    bool stopping = false;

public:
    explicit ThreadPool(size_t threadCount)
    {
        if (threadCount == 0)
        {
            throw std::invalid_argument("Thread count must be positive");
        }

        for (size_t i = 0; i < threadCount; ++i)
        {

            workers.emplace_back([this]()
                                 {

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
                } });
        }
    }

    template <typename F>
    auto submit(F &&function)
        -> std::future<std::invoke_result_t<F>>
    {

        using ReturnType = std::invoke_result_t<F>;

        auto task = std::make_shared<
            std::packaged_task<ReturnType()>>(std::forward<F>(function));

        std::future<ReturnType> result = task->get_future();

        {
            std::lock_guard<std::mutex> lock(queueMutex);

            if (stopping)
            {
                throw std::runtime_error("Thread pool is stopping");
            }

            tasks.emplace([task]()
                          { (*task)(); });
        }

        condition.notify_one();

        return result;
    }

    ~ThreadPool()
    {

        {
            std::lock_guard<std::mutex> lock(queueMutex);
            stopping = true;
        }

        condition.notify_all();

        for (auto &worker : workers)
        {
            if (worker.joinable())
            {
                worker.join();
            }
        }
    }

    ThreadPool(const ThreadPool &) = delete;
    ThreadPool &operator=(const ThreadPool &) = delete;
};

// High-quality hash-based PRNG returning uniform double in [0, 1)
// Function of (seed, day, cell, stream) as specified in Section 3.1
inline uint64_t hash_coords(uint64_t seed, uint64_t day, uint64_t cell, uint64_t stream)
{
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

inline double rand01(uint64_t seed, uint64_t day, uint64_t cell, uint64_t stream)
{
    uint64_t h = hash_coords(seed, day, cell, stream);
    return (h >> 11) * (1.0 / 9007199254740992.0); // 1.0 / 2^53
}

// 64-bit FNV-1a hash over the entire grid
uint64_t computeGridChecksum(const std::vector<Cell> &grid)
{
    uint64_t hash = 14695981039346656037ULL;
    for (const auto &cell : grid)
    {
        uint16_t val = (static_cast<uint16_t>(cell.state) << 8) | cell.c;
        hash ^= (val & 0xFF);
        hash *= 1099511628211ULL;
        hash ^= ((val >> 8) & 0xFF);
        hash *= 1099511628211ULL;
    }
    return hash;
}

Counts updateRegion(
    const std::vector<Cell> &cur,
    std::vector<Cell> &next,
    const Region &region,
    int N,
    int day,
    const double P[9],
    const Params &p)
{
    Counts localCounts;

    const int di[8] = {-1, -1, -1, 0, 0, 1, 1, 1};
    const int dj[8] = {-1, 0, 1, -1, 1, -1, 0, 1};

    for (int i = region.rowStart;
         i < region.rowEnd;
         i += region.rowStride)
    {
        for (int j = region.colStart; j < region.colEnd; ++j)
        {

            long long id = 1LL * i * N + j;
            const Cell &cell = cur[id];

            if (cell.state == STATE_S)
            {

                int k = 0;

                for (int dir = 0; dir < 8; ++dir)
                {
                    int ni = i + di[dir];
                    int nj = j + dj[dir];

                    if (ni >= 0 && ni < N &&
                        nj >= 0 && nj < N &&
                        cur[1LL * ni * N + nj].state == STATE_I)
                    {
                        k++;
                    }
                }

                if (k > 0 &&
                    rand01(p.seed, day, id, 0) < P[k])
                {
                    next[id] = {STATE_I, 0};
                    localCounts.I++;
                }
                else
                {
                    next[id] = {STATE_S, 0};
                    localCounts.S++;
                }
            }
            else if (cell.state == STATE_I)
            {

                if (cell.c + 1 == p.d)
                {

                    if (rand01(p.seed, day, id, 1) < p.mu)
                    {
                        next[id] = {STATE_D, 0};
                        localCounts.D++;
                    }
                    else
                    {
                        next[id] = {STATE_R, 0};
                        localCounts.R++;
                    }
                }
                else
                {
                    next[id] = {
                        STATE_I,
                        static_cast<uint8_t>(cell.c + 1)};
                    localCounts.I++;
                }
            }
            else
            {

                next[id] = cell;

                if (cell.state == STATE_R)
                    localCounts.R++;
                else if (cell.state == STATE_V)
                    localCounts.V++;
                else if (cell.state == STATE_D)
                    localCounts.D++;
            }
        }
    }

    return localCounts;
}

// Keep this wrapper for existing B1 and B2 code
Counts updateChunk(
    const std::vector<Cell> &cur,
    std::vector<Cell> &next,
    int startRow,
    int endRow,
    int N,
    int day,
    const double P[9],
    const Params &p)
{
    Region region{startRow, endRow, 0, N};

    return updateRegion(
        cur, next, region, N, day, P, p);
}

std::string snapshotFilename(const std::string &prefix, int day)
{
    std::ostringstream filename;
    filename << prefix << "_day_"
             << std::setw(4) << std::setfill('0')
             << day << ".ppm";
    return filename.str();
}

bool writePPM(
    const std::vector<Cell> &grid,
    int N,
    const std::string &filename,
    int latencyMs)
{
    try
    {
        std::filesystem::path outputPath(filename);

        if (!outputPath.parent_path().empty())
        {
            std::filesystem::create_directories(
                outputPath.parent_path());
        }

        if (latencyMs > 0)
        {
            std::this_thread::sleep_for(
                std::chrono::milliseconds(latencyMs));
        }

        std::ofstream out(filename, std::ios::binary);

        if (!out.is_open())
        {
            return false;
        }

        out << "P6\n"
            << N << " " << N << "\n255\n";

        for (const auto &cell : grid)
        {
            unsigned char rgb[3];

            switch (cell.state)
            {
            case STATE_S:
                rgb[0] = 245;
                rgb[1] = 245;
                rgb[2] = 245;
                break;

            case STATE_I:
                rgb[0] = 220;
                rgb[1] = 40;
                rgb[2] = 40;
                break;

            case STATE_R:
                rgb[0] = 40;
                rgb[1] = 180;
                rgb[2] = 80;
                break;

            case STATE_V:
                rgb[0] = 40;
                rgb[1] = 100;
                rgb[2] = 220;
                break;

            case STATE_D:
            default:
                rgb[0] = 0;
                rgb[1] = 0;
                rgb[2] = 0;
                break;
            }

            out.write(
                reinterpret_cast<const char *>(rgb), 3);
        }

        out.flush();
        bool success = out.good();
        out.close();

        return success && !out.fail();
    }
    catch (...)
    {
        return false;
    }
}

SweepSummary runSweepSimulation(
    const Params &base,
    const SweepSpec &spec)
{
    Params p = base;

    p.beta = spec.beta;
    p.v = spec.vaccination;
    p.seed = spec.seed;

    const int N = p.N;
    const long long total_cells = 1LL * N * N;

    std::vector<Cell> cur(total_cells, {STATE_S, 0});
    std::vector<Cell> next(total_cells, {STATE_S, 0});

    std::mt19937_64 rng(p.seed);

    // Vaccination initialization
    long long num_vaccinated =
        static_cast<long long>(
            std::round(p.v * total_cells));

    num_vaccinated = std::clamp(
        num_vaccinated, 0LL, total_cells);

    if (num_vaccinated > 0)
    {
        std::vector<long long> indices(total_cells);

        std::iota(indices.begin(), indices.end(), 0LL);

        for (long long i = 0; i < num_vaccinated; ++i)
        {
            std::uniform_int_distribution<long long>
                dist(i, total_cells - 1);

            long long j = dist(rng);

            std::swap(indices[i], indices[j]);

            cur[indices[i]].state = STATE_V;
        }
    }

    // Initial infection
    long long susceptible =
        total_cells - num_vaccinated;

    long long initial_infected =
        std::min(static_cast<long long>(p.I0), susceptible);

    long long infected_so_far = 0;

    while (infected_so_far < initial_infected)
    {
        std::uniform_int_distribution<long long>
            dist(0, total_cells - 1);

        long long idx = dist(rng);

        if (cur[idx].state == STATE_S)
        {
            cur[idx].state = STATE_I;
            cur[idx].c = 0;
            infected_so_far++;
        }
    }

    // Initial counts
    Counts counts;

    for (const auto &cell : cur)
    {
        switch (cell.state)
        {
        case STATE_S:
            counts.S++;
            break;
        case STATE_I:
            counts.I++;
            break;
        case STATE_R:
            counts.R++;
            break;
        case STATE_V:
            counts.V++;
            break;
        case STATE_D:
            counts.D++;
            break;
        }
    }

    long long total_infected = initial_infected;
    long long peak_infected = counts.I;

    int peak_day = 0;
    int days_run = 0;

    bool lockdown = false;

    const Region wholeGrid{0, N, 0, N, 1};

    for (int day = 0; day < p.max_days; ++day)
    {

        if (!p.fixed_days && counts.I == 0)
        {
            break;
        }

        double infected_fraction =
            static_cast<double>(counts.I) /
            static_cast<double>(total_cells);

        if (!lockdown)
        {
            if (infected_fraction > p.theta_on)
            {
                lockdown = true;
            }
        }
        else
        {
            if (infected_fraction < p.theta_off)
            {
                lockdown = false;
            }
        }

        double beta_eff =
            p.beta * p.m *
            (lockdown ? p.lambda : 1.0);

        double P[9];

        for (int k = 0; k <= 8; ++k)
        {
            P[k] = 1.0 - std::pow(1.0 - beta_eff, k);
        }

        // Sequential execution inside this independent run
        Counts next_counts = updateRegion(
            cur, next, wholeGrid,
            N, day, P, p);

        // Count newly infected people
        for (long long i = 0; i < total_cells; ++i)
        {
            if (cur[i].state == STATE_S &&
                next[i].state == STATE_I)
            {
                total_infected++;
            }
        }

        std::swap(cur, next);

        counts = next_counts;
        days_run++;

        if (counts.I > peak_infected)
        {
            peak_infected = counts.I;
            peak_day = day + 1;
        }
    }

    SweepSummary result;

    result.run_id = spec.run_id;
    result.beta = spec.beta;
    result.vaccination = spec.vaccination;
    result.seed = spec.seed;

    result.peak_infected_fraction =
        static_cast<double>(peak_infected) /
        static_cast<double>(total_cells);

    result.day_of_peak = peak_day;
    result.total_infected = total_infected;
    result.total_deaths = counts.D;
    result.days = days_run;

    return result;
}

int runParameterSweep(const Params &p)
{

    if (p.N <= 0 || p.max_days < 0 ||
        p.sweep_workers <= 0 ||
        p.time_limit_seconds < 0 ||
        p.sweep_betas.empty() ||
        p.sweep_vaccinations.empty() ||
        p.sweep_seeds.empty())
    {

        std::cerr << "Invalid parameter sweep settings\n";
        return 1;
    }

    // Generate all independent parameter combinations
    std::vector<SweepSpec> specs;

    int runId = 1;

    for (double beta : p.sweep_betas)
    {
        for (double vaccination : p.sweep_vaccinations)
        {
            for (uint64_t seed : p.sweep_seeds)
            {

                specs.push_back({runId++, beta, vaccination, seed});
            }
        }
    }

    const size_t totalRuns = specs.size();

    std::filesystem::path outputPath(p.sweep_output);

    if (!outputPath.parent_path().empty())
    {
        std::filesystem::create_directories(
            outputPath.parent_path());
    }

    std::ofstream csv(p.sweep_output);

    if (!csv.is_open())
    {
        std::cerr << "Cannot open sweep output file\n";
        return 1;
    }

    csv << "run_id,beta,vaccination,seed,"
        << "peak_infected_fraction,day_of_peak,"
        << "total_infected,total_deaths,days\n";

    struct ActiveRun
    {
        size_t index;
        std::future<SweepSummary> future;
    };

    std::vector<ActiveRun> active;

    std::queue<size_t> completedQueue;
    std::mutex completedMutex;
    std::condition_variable completedCondition;

    size_t nextIndex = 0;
    size_t finished = 0;
    size_t failed = 0;

    const auto startTime = std::chrono::steady_clock::now();

    auto timeExpired = [&]()
    {
        if (p.time_limit_seconds == 0)
        {
            return false;
        }

        return std::chrono::steady_clock::now() - startTime >= std::chrono::seconds(p.time_limit_seconds);
    };

    std::cout << "Total independent runs: "
              << totalRuns << "\n";

    while (nextIndex < totalRuns || !active.empty())
    {

        // Submit new runs while workers are available
        while (nextIndex < totalRuns &&
               active.size() < static_cast<size_t>(p.sweep_workers) &&
               !timeExpired())
        {

            size_t index = nextIndex++;
            SweepSpec spec = specs[index];

            auto future = std::async(
                std::launch::async,
                [&, index, spec, p]()
                {
                    try
                    {
                        SweepSummary summary =
                            runSweepSimulation(p, spec);

                        {
                            std::lock_guard<std::mutex> lock(
                                completedMutex);

                            completedQueue.push(index);
                        }

                        completedCondition.notify_one();

                        return summary;
                    }
                    catch (...)
                    {

                        {
                            std::lock_guard<std::mutex> lock(
                                completedMutex);

                            completedQueue.push(index);
                        }

                        completedCondition.notify_one();

                        throw;
                    }
                });

            active.push_back({index, std::move(future)});

            std::cout << "Started run "
                      << spec.run_id << "\n";
        }

        if (active.empty())
        {
            break;
        }

        // Receive whichever run finishes first
        size_t completedIndex;

        {
            std::unique_lock<std::mutex> lock(
                completedMutex);

            completedCondition.wait(
                lock,
                [&]()
                {
                    return !completedQueue.empty();
                });

            completedIndex = completedQueue.front();
            completedQueue.pop();
        }

        auto it = std::find_if(
            active.begin(),
            active.end(),
            [completedIndex](const ActiveRun &run)
            {
                return run.index == completedIndex;
            });

        if (it == active.end())
        {
            std::cerr << "Internal sweep scheduling error\n";
            return 1;
        }

        try
        {

            SweepSummary result = it->future.get();

            csv << result.run_id << ","
                << std::fixed << std::setprecision(6)
                << result.beta << ","
                << result.vaccination << ","
                << result.seed << ","
                << result.peak_infected_fraction << ","
                << result.day_of_peak << ","
                << result.total_infected << ","
                << result.total_deaths << ","
                << result.days << "\n";

            csv.flush();
        }
        catch (const std::exception &e)
        {

            failed++;

            std::cerr << "Run "
                      << specs[completedIndex].run_id
                      << " failed: " << e.what() << "\n";
        }

        active.erase(it);

        finished++;

        std::cout << finished << "/"
                  << totalRuns << " runs done\n";
    }

    // Report runs that were never started
    if (nextIndex < totalRuns)
    {

        std::cout << "\nMissing runs (time limit reached):\n";

        for (size_t i = nextIndex; i < totalRuns; ++i)
        {

            std::cout << "Run " << specs[i].run_id
                      << " | beta=" << specs[i].beta
                      << " | vaccination=" << specs[i].vaccination
                      << " | seed=" << specs[i].seed
                      << "\n";
        }
    }

    csv.close();

    std::cout << "\nParameter sweep completed.\n";
    std::cout << "Finished: " << finished << "/"
              << totalRuns << "\n";
    std::cout << "Failed: " << failed << "\n";
    std::cout << "Not started: "
              << totalRuns - nextIndex << "\n";
    std::cout << "Results: " << p.sweep_output << "\n";

    return failed == 0 ? 0 : 1;
}

void printHelp(const char *prog)
{
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
              << "  --async-snapshots             Enable asynchronous snapshots\n"
              << "  --sweep                       Run parameter sweep\n"
              << "  --sweep-betas <list>          Comma-separated beta values\n"
              << "  --sweep-vaccinations <list>   Comma-separated vaccination values\n"
              << "  --sweep-seeds <list>          Comma-separated seeds\n"
              << "  --sweep-workers <int>         Maximum concurrent runs\n"
              << "  --time-limit <seconds>        Sweep time limit; 0 = unlimited\n"
              << "  --sweep-output <file>         Sweep results CSV\n"
              << "  --snapshot-prefix <path>      Snapshot filename prefix\n"
              << "  --fixed-days                  Simulate exact number of days even if I reaches 0\n"
              << "  -o, --output <file>           Output CSV filename (default: output.csv)\n"
              << "  -h, --help                    Show this help message\n"
              << "  --pool                       Use reusable thread pool\n"
              << "  --partition <type>    strips, blocks or interleaved\n"
              << "  --sweep                       Run parameter sweep\n"
              << "  --sweep-betas <list>          Beta values\n"
              << "  --sweep-vaccinations <list>   Vaccination values\n"
              << "  --sweep-seeds <list>          Random seeds\n"
              << "  --sweep-workers <int>         Concurrent runs\n"
              << "  --time-limit <seconds>       Time limit\n"
              << "  --sweep-output <file>         Results CSV\n";
}

std::vector<double> parseDoubleList(const std::string &input)
{

    std::vector<double> values;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
    {

        if (token.empty())
        {
            throw std::invalid_argument("Invalid numeric list");
        }

        double value = std::stod(token);

        if (!std::isfinite(value) || value < 0.0 || value > 1.0)
        {
            throw std::invalid_argument("Values must be between 0 and 1");
        }

        values.push_back(value);
    }

    if (values.empty())
    {
        throw std::invalid_argument("Empty numeric list");
    }

    return values;
}

std::vector<uint64_t> parseSeedList(const std::string &input)
{

    std::vector<uint64_t> values;
    std::stringstream ss(input);
    std::string token;

    while (std::getline(ss, token, ','))
    {

        if (token.empty() || token[0] == '-')
        {
            throw std::invalid_argument("Invalid seed list");
        }

        values.push_back(std::stoull(token));
    }

    if (values.empty())
    {
        throw std::invalid_argument("Empty seed list");
    }

    return values;
}

Params parseArgs(int argc, char *argv[])
{
    Params p;
    for (int i = 1; i < argc; ++i)
    {
        std::string arg = argv[i];
        if ((arg == "-N" || arg == "--grid-size") && i + 1 < argc)
        {
            p.N = std::stoi(argv[++i]);
        }
        else if ((arg == "-D" || arg == "--days") && i + 1 < argc)
        {
            p.max_days = std::stoi(argv[++i]);
        }
        else if ((arg == "-s" || arg == "--seed") && i + 1 < argc)
        {
            p.seed = std::stoull(argv[++i]);
        }
        else if ((arg == "-b" || arg == "--beta") && i + 1 < argc)
        {
            p.beta = std::stod(argv[++i]);
        }
        else if ((arg == "-d" || arg == "--days-infectious") && i + 1 < argc)
        {
            p.d = std::stoi(argv[++i]);
        }
        else if ((arg == "-u" || arg == "--mu" || arg == "--mortality") && i + 1 < argc)
        {
            p.mu = std::stod(argv[++i]);
        }
        else if ((arg == "--I0" || arg == "--initially-infected") && i + 1 < argc)
        {
            p.I0 = std::stoi(argv[++i]);
        }
        else if ((arg == "-v" || arg == "--vaccinated") && i + 1 < argc)
        {
            p.v = std::stod(argv[++i]);
        }
        else if ((arg == "-m" || arg == "--mask-factor") && i + 1 < argc)
        {
            p.m = std::stod(argv[++i]);
        }
        else if ((arg == "-l" || arg == "--lambda") && i + 1 < argc)
        {
            p.lambda = std::stod(argv[++i]);
        }
        else if (arg == "--theta-on" && i + 1 < argc)
        {
            p.theta_on = std::stod(argv[++i]);
        }
        else if (arg == "--theta-off" && i + 1 < argc)
        {
            p.theta_off = std::stod(argv[++i]);
        }
        else if ((arg == "-T" || arg == "--threads") && i + 1 < argc)
        {
            p.num_threads = std::stoi(argv[++i]);
        }
        else if ((arg == "-C" || arg == "--chunks") && i + 1 < argc)
        {
            p.num_chunks = std::stoi(argv[++i]);
        }
        else if ((arg == "-K" || arg == "--snapshot-interval") && i + 1 < argc)
        {
            p.snapshot_interval = std::stoi(argv[++i]);
        }
        else if ((arg == "-Q" || arg == "--snapshot-limit") && i + 1 < argc)
        {
            p.snapshot_limit = std::stoi(argv[++i]);
        }
        else if (arg == "--io-latency" && i + 1 < argc)
        {
            p.io_latency_ms = std::stoi(argv[++i]);
        }
        else if (arg == "--async-snapshots")
        {
            p.async_snapshots = true;
        }
        else if (arg == "--snapshot-prefix" && i + 1 < argc)
        {
            p.snapshot_prefix = argv[++i];
        }
        else if (arg == "--sweep")
        {
            p.run_sweep = true;
        }
        else if (arg == "--sweep-betas" && i + 1 < argc)
        {
            p.sweep_betas = parseDoubleList(argv[++i]);
        }
        else if (arg == "--sweep-vaccinations" && i + 1 < argc)
        {
            p.sweep_vaccinations = parseDoubleList(argv[++i]);
        }
        else if (arg == "--sweep-seeds" && i + 1 < argc)
        {
            p.sweep_seeds = parseSeedList(argv[++i]);
        }
        else if (arg == "--sweep-workers" && i + 1 < argc)
        {
            p.sweep_workers = std::stoi(argv[++i]);
        }
        else if (arg == "--time-limit" && i + 1 < argc)
        {
            p.time_limit_seconds = std::stoi(argv[++i]);
        }
        else if (arg == "--sweep-output" && i + 1 < argc)
        {
            p.sweep_output = argv[++i];
        }
        else if (arg == "--sweep")
        {
            p.run_sweep = true;
        }
        else if (arg == "--sweep-betas" && i + 1 < argc)
        {
            p.sweep_betas = parseDoubleList(argv[++i]);
        }
        else if (arg == "--sweep-vaccinations" && i + 1 < argc)
        {
            p.sweep_vaccinations = parseDoubleList(argv[++i]);
        }
        else if (arg == "--sweep-seeds" && i + 1 < argc)
        {
            p.sweep_seeds = parseSeedList(argv[++i]);
        }
        else if (arg == "--sweep-workers" && i + 1 < argc)
        {
            p.sweep_workers = std::stoi(argv[++i]);
        }
        else if (arg == "--time-limit" && i + 1 < argc)
        {
            p.time_limit_seconds = std::stoi(argv[++i]);
        }
        else if (arg == "--sweep-output" && i + 1 < argc)
        {
            p.sweep_output = argv[++i];
        }
        else if (arg == "--fixed-days")
        {
            p.fixed_days = true;
        }
        else if ((arg == "-o" || arg == "--output") && i + 1 < argc)
        {
            p.output_file = argv[++i];
        }
        else if (arg == "-h" || arg == "--help")
        {
            printHelp(argv[0]);
            std::exit(0);
        }
        else if (arg == "--pool")
        {
            p.use_thread_pool = true;
        }
        else if (arg == "--partition" && i + 1 < argc)
        {
            p.partition = argv[++i];
        }
        else
        {
            std::cerr << "Warning: Unknown option " << arg << "\n";
        }
    }
    return p;
}

int main(int argc, char *argv[])
{
    Params p = parseArgs(argc, argv);
    if (p.run_sweep)
    {
        return runParameterSweep(p);
    }
    if (p.snapshot_interval < 0 || p.io_latency_ms < 0 ||
        (p.async_snapshots && p.snapshot_interval > 0 &&
         p.snapshot_limit <= 0))
    {

        std::cerr << "Invalid snapshot parameters\n";
        return 1;
    }

    const int N = p.N;
    const long long total_cells = static_cast<long long>(N) * N;

    // Grid buffers: current and next
    std::vector<Cell> cur(total_cells, {STATE_S, 0});
    std::vector<Cell> next(total_cells, {STATE_S, 0});

    // Initial conditions setup (done sequentially)
    std::mt19937_64 rng(p.seed);

    // 1. Vaccinate fraction v
    long long num_vaccinated = static_cast<long long>(std::round(p.v * total_cells));
    if (num_vaccinated > total_cells)
        num_vaccinated = total_cells;

    if (num_vaccinated > 0)
    {
        std::vector<long long> all_indices(total_cells);
        std::iota(all_indices.begin(), all_indices.end(), 0LL);
        // Partial Fisher-Yates shuffle
        for (long long i = 0; i < num_vaccinated; ++i)
        {
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

    while (infected_so_far < to_infect)
    {
        std::uniform_int_distribution<long long> dist(0, total_cells - 1);
        long long idx = dist(rng);
        if (cur[idx].state == STATE_S)
        {
            cur[idx].state = STATE_I;
            cur[idx].c = 0;
            infected_so_far++;
        }
    }

    // Open output CSV
    std::ofstream csv(p.output_file);
    if (!csv.is_open())
    {
        std::cerr << "Error: Could not open output file " << p.output_file << "\n";
        return 1;
    }
    csv << "day,S,I,R,V,D,beta_eff,checksum\n";

    // Initial counts for day 0
    Counts counts;
    for (long long idx = 0; idx < total_cells; ++idx)
    {
        switch (cur[idx].state)
        {
        case STATE_S:
            counts.S++;
            break;
        case STATE_I:
            counts.I++;
            break;
        case STATE_R:
            counts.R++;
            break;
        case STATE_V:
            counts.V++;
            break;
        case STATE_D:
            counts.D++;
            break;
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
    std::deque<std::pair<int, std::future<bool>>> pendingSnapshots;

    if (p.use_thread_pool)
    {
        pool = std::make_unique<ThreadPool>(p.num_threads);
    }

    // Main simulation loop
    for (int day = 0; day < p.max_days; ++day)
    {
        // Check stopping condition: stop if no infected and not fixed_days
        if (!p.fixed_days && counts.I == 0)
        {
            break;
        }

        // 1. Determine beta_eff for the upcoming day
        // Decisions made from day t's totals apply to day t + 1
        double infected_fraction = static_cast<double>(counts.I) / static_cast<double>(total_cells);
        if (!lockdown)
        {
            if (infected_fraction > p.theta_on)
            {
                lockdown = true;
            }
        }
        else
        {
            if (infected_fraction < p.theta_off)
            {
                lockdown = false;
            }
        }

        beta_eff = (p.beta * p.m) * (lockdown ? p.lambda : 1.0);

        // Precompute P[k] = 1 - (1 - beta_eff)^k for k = 0..8
        double P[9];
        for (int k = 0; k <= 8; ++k)
        {
            P[k] = 1.0 - std::pow(1.0 - beta_eff, k);
        }

        // Reset next counts
        Counts next_counts;

        int C = std::max(1, p.num_chunks);
        C = std::min(C, N);

        std::vector<Region> regions;
        regions.reserve(C);

        if (p.partition == "strips")
        {

            // B3.1: Row strips
            int rowsPerChunk = N / C;
            int extraRows = N % C;
            int startRow = 0;

            for (int chunk = 0; chunk < C; ++chunk)
            {
                int rows = rowsPerChunk +
                           (chunk < extraRows ? 1 : 0);

                int endRow = startRow + rows;

                regions.push_back({startRow, endRow, 0, N});

                startRow = endRow;
            }
        }
        else if (p.partition == "blocks")
        {

            // B3.2: Approximately square blocks
            int tileRows = static_cast<int>(std::ceil(std::sqrt(C)));

            int baseBlocks = C / tileRows;
            int extraBlocks = C % tileRows;

            int rowStart = 0;

            for (int r = 0; r < tileRows; ++r)
            {

                int rows = N / tileRows +
                           (r < N % tileRows ? 1 : 0);

                int rowEnd = rowStart + rows;

                int blocksInRow = baseBlocks +
                                  (r < extraBlocks ? 1 : 0);

                int colStart = 0;

                for (int c = 0; c < blocksInRow; ++c)
                {

                    int cols = N / blocksInRow +
                               (c < N % blocksInRow ? 1 : 0);

                    int colEnd = colStart + cols;

                    regions.push_back({rowStart, rowEnd,
                                       colStart, colEnd});

                    colStart = colEnd;
                }

                rowStart = rowEnd;
            }
        }
        else if (p.partition == "interleaved")
        {

            // B3.3: Interleaved rows
            int T = p.num_threads;

            if (T <= 0 || C < T)
            {
                throw std::invalid_argument(
                    "Interleaved partition requires C >= T > 0");
            }

            for (int chunk = 0; chunk < C; ++chunk)
            {

                int lane = chunk % T;
                int subtask = chunk / T;

                int tasksForLane = 1 + (C - 1 - lane) / T;

                int rowStart = lane + subtask * T;
                int stride = T * tasksForLane;

                regions.push_back({rowStart,
                                   N,
                                   0,
                                   N,
                                   stride});
            }
        }
        else
        {
            throw std::invalid_argument(
                "Unknown partition. Use strips, blocks or interleaved.");
        }

        // Execute generated regions
        next_counts = Counts{};

        if (pool)
        {

            // B2: Thread pool
            std::vector<std::future<Counts>> futures;
            futures.reserve(regions.size());

            for (size_t chunk = 0; chunk < regions.size(); ++chunk)
            {

                futures.emplace_back(
                    pool->submit(
                        [&, chunk]()
                        {
                            return updateRegion(
                                cur, next, regions[chunk],
                                N, day, P, p);
                        }));
            }

            for (auto &future : futures)
            {
                Counts local = future.get();

                next_counts.S += local.S;
                next_counts.I += local.I;
                next_counts.R += local.R;
                next_counts.V += local.V;
                next_counts.D += local.D;
            }
        }
        else
        {

            // B1: Thread per day
            std::vector<std::thread> workers;
            std::vector<Counts> threadCounts(regions.size());

            for (size_t chunk = 0; chunk < regions.size(); ++chunk)
            {

                workers.emplace_back(
                    [&, chunk]()
                    {
                        threadCounts[chunk] = updateRegion(
                            cur, next, regions[chunk],
                            N, day, P, p);
                    });
            }

            for (auto &worker : workers)
            {
                worker.join();
            }

            for (const auto &local : threadCounts)
            {
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

        // C1: Background snapshot writer

        if (p.snapshot_interval > 0 &&
            (day + 1) % p.snapshot_interval == 0)
        {

            int snapshotDay = day + 1;

            std::string filename = snapshotFilename(
                p.snapshot_prefix, snapshotDay);

            if (!p.async_snapshots)
            {

                // Blocking snapshot
                if (!writePPM(
                        cur, N, filename, p.io_latency_ms))
                {
                    std::cerr << "Snapshot failed: "
                              << filename << "\n";
                }
            }
            else
            {

                // Wait for oldest write if Q writes are pending
                if (pendingSnapshots.size() >=
                    static_cast<size_t>(p.snapshot_limit))
                {

                    auto oldest = std::move(
                        pendingSnapshots.front());

                    pendingSnapshots.pop_front();

                    if (!oldest.second.get())
                    {
                        std::cerr << "Snapshot failed on day "
                                  << oldest.first << "\n";
                    }
                }

                // Copy current grid before submitting async task
                std::vector<Cell> snapshot = cur;

                auto future = std::async(
                    std::launch::async,
                    [snapshot = std::move(snapshot),
                     N, filename,
                     latency = p.io_latency_ms]()
                    {
                        return writePPM(
                            snapshot, N, filename, latency);
                    });

                pendingSnapshots.emplace_back(
                    snapshotDay, std::move(future));
            }
        }
    }
    // Wait for all remaining asynchronous snapshot writes
    while (!pendingSnapshots.empty())
    {

        auto task = std::move(pendingSnapshots.front());
        pendingSnapshots.pop_front();

        if (!task.second.get())
        {
            std::cerr << "Snapshot failed on day "
                      << task.first << "\n";
        }
    }

    csv.close();
    std::cout << "Simulation completed. Results written to " << p.output_file << "\n";
    return 0;
}
