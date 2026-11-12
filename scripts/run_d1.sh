# ARYAN SHOKEEN MT26116
#!/bin/bash
set -euo pipefail

mkdir -p results/part_d

make

LOGICAL_CORES=$(sysctl -n hw.logicalcpu)
MAX_THREADS=$((2 * LOGICAL_CORES))

THREADS=(1)
T=2

while [ "$T" -le "$MAX_THREADS" ]; do
    THREADS+=("$T")
    T=$((T * 2))
done

RAW="results/part_d/d1_raw.csv"

echo "threads,run,seconds" > "$RAW"

echo "Logical cores: $LOGICAL_CORES"
echo "Maximum threads: $MAX_THREADS"
echo "Testing threads: ${THREADS[*]}"

for T in "${THREADS[@]}"; do

    for RUN in 1 2 3 4 5; do

        echo "Running T=$T, repetition=$RUN"

        OUTPUT=$(./pandemic_sim \
            -N 4000 \
            -D 200 \
            -s 42 \
            -T "$T" \
            -C "$T" \
            --partition strips \
            --fixed-days \
            --benchmark)

        echo "$OUTPUT"

        ELAPSED_TIME=$(echo "$OUTPUT" | awk '/Simulation loop time:/ {print $4}')

        if [ -z "$SECONDS" ]; then
            echo "Could not retrieve execution time"
            exit 1
        fi

        echo "$T,$RUN,$ELAPSED_TIME" >> "$RAW"

    done
done

echo "D1 benchmark completed."
echo "Raw results: $RAW"

