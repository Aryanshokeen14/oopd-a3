
#!/bin/bash
set -euo pipefail

mkdir -p results/part_d
make

T=8
N=4000
DAYS=200
REPEATS=5

CHUNKS=(8 16 32 64 128 256 512)

RAW="results/part_d/d3_chunks_raw.csv"

echo "threads,chunks,run,seconds" > "$RAW"

for C in "${CHUNKS[@]}"; do

    for RUN in $(seq 1 "$REPEATS"); do

        echo "Testing T=$T C=$C repetition=$RUN"

        OUTPUT=$(./pandemic_sim \
            -N "$N" -D "$DAYS" -s 42 \
            -T "$T" -C "$C" \
            --pool --partition strips \
            --fixed-days --benchmark)

        ELAPSED_TIME=$(echo "$OUTPUT" |
            awk '/Simulation loop time:/ {print $4}')

        if [ -z "$ELAPSED_TIME" ] || [ "$ELAPSED_TIME" = "0" ]; then
            echo "Invalid timing"
            exit 1
        fi

        echo "$OUTPUT"

        echo "$T,$C,$RUN,$ELAPSED_TIME" >> "$RAW"

    done
done

echo "D3.1 benchmark completed."

