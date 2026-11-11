
#!/bin/bash
set -euo pipefail

mkdir -p results/part_d
make

THREADS=8
CHUNKS=8
DAYS=200
REPEATS=5

SIZES=(50 100 200 500 1000)
MODES=(sequential thread-per-day thread-pool)

RAW="results/part_d/d3_small_grids_raw.csv"

echo "N,mode,run,seconds" > "$RAW"

for N in "${SIZES[@]}"; do

    for MODE in "${MODES[@]}"; do

        for RUN in $(seq 1 "$REPEATS"); do

            echo "Testing N=$N Mode=$MODE Run=$RUN"

            if [ "$MODE" = "sequential" ]; then

                OUTPUT=$(./pandemic_sim \
                    -N "$N" -D "$DAYS" -s 42 \
                    --fixed-days --benchmark --sequential)

            elif [ "$MODE" = "thread-per-day" ]; then

                OUTPUT=$(./pandemic_sim \
                    -N "$N" -D "$DAYS" -s 42 \
                    -T "$THREADS" -C "$CHUNKS" \
                    --partition strips \
                    --fixed-days --benchmark)

            else

                OUTPUT=$(./pandemic_sim \
                    -N "$N" -D "$DAYS" -s 42 \
                    -T "$THREADS" -C "$CHUNKS" \
                    --partition strips \
                    --pool --fixed-days --benchmark)

            fi

            ELAPSED_TIME=$(echo "$OUTPUT" |
                awk '/Simulation loop time:/ {print $4}')

            if [ -z "$ELAPSED_TIME" ] || [ "$ELAPSED_TIME" = "0" ]; then
                echo "Invalid timing"
                exit 1
            fi

            echo "$OUTPUT"

            echo "$N,$MODE,$RUN,$ELAPSED_TIME" >> "$RAW"

        done
    done
done

echo "D3.3 benchmark completed."

