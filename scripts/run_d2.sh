
#!/bin/bash
set -euo pipefail

mkdir -p results/part_d
make

THREADS=(1 2 4 8 16)
SIZES=(1000 1414 2000 2828 4000)

RAW="results/part_d/d2_raw.csv"

echo "threads,grid_size,run,seconds" > "$RAW"

for INDEX in "${!THREADS[@]}"; do

    T=${THREADS[$INDEX]}
    N=${SIZES[$INDEX]}

    for RUN in 1 2 3 4 5; do

        echo "Running T=$T N=$N repetition=$RUN"

        OUTPUT=$(./pandemic_sim \
            -N "$N" \
            -D 200 \
            -s 42 \
            -T "$T" \
            -C "$T" \
            --partition strips \
            --fixed-days \
            --benchmark)

        echo "$OUTPUT"

        ELAPSED_TIME=$(echo "$OUTPUT" |
            awk '/Simulation loop time:/ {print $4}')

        if [ -z "$ELAPSED_TIME" ]; then
            echo "Could not retrieve execution time"
            exit 1
        fi

        echo "$T,$N,$RUN,$ELAPSED_TIME" >> "$RAW"

    done
done

echo "D2 benchmark completed."
echo "Results: $RAW"

