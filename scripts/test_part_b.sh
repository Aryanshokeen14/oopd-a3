# ARYAN SHOKEEN MT26116
#!/bin/bash
set -euo pipefail

mkdir -p results/part_b

make

TMP_DIR=$(mktemp -d)
trap 'rm -rf "$TMP_DIR"' EXIT

# Compile original sequential Part A
git show main:pandemic_sim.cpp > "$TMP_DIR/pandemic_A.cpp"

clang++ -O3 -std=c++17 -Wall -Wextra \
    "$TMP_DIR/pandemic_A.cpp" \
    -o "$TMP_DIR/pandemic_A"

# Generate sequential reference
"$TMP_DIR/pandemic_A" \
    -N 100 -D 50 -s 42 \
    -o results/part_b/sequential.csv \
    --fixed-days

# Required thread counts
THREADS=(1 2 3 4 7 8 16)

# All partitioning strategies
PARTITIONS=(strips blocks interleaved)

for T in "${THREADS[@]}"; do

    for PARTITION in "${PARTITIONS[@]}"; do

        for MODE in thread-per-day thread-pool; do

            if [ "$MODE" = "thread-pool" ]; then
                C=$((T * 2))
                POOL_OPTION="--pool"
            else
                C=$T
                POOL_OPTION=""
            fi

            OUTPUT="results/part_b/${MODE}_T${T}_C${C}_${PARTITION}.csv"

            echo "Testing: T=$T C=$C Mode=$MODE Partition=$PARTITION"

            if [ "$MODE" = "thread-pool" ]; then
                ./pandemic_sim \
                    -N 100 -D 50 -s 42 \
                    -T "$T" -C "$C" \
                    --partition "$PARTITION" \
                    --fixed-days --pool \
                    -o "$OUTPUT"
            else
                ./pandemic_sim \
                    -N 100 -D 50 -s 42 \
                    -T "$T" -C "$C" \
                    --partition "$PARTITION" \
                    --fixed-days \
                    -o "$OUTPUT"
            fi

            if diff -q results/part_b/sequential.csv "$OUTPUT" >/dev/null; then
                echo "PASS"
            else
                echo "FAIL: $OUTPUT"
                diff -u results/part_b/sequential.csv "$OUTPUT"
                exit 1
            fi

        done
    done
done

echo "=================================="
echo "ALL PART B CORRECTNESS TESTS PASSED"
echo "=================================="

