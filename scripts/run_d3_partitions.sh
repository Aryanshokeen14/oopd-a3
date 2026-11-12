# ARYAN SHOKEEN MT26116
#!/bin/bash
set -euo pipefail

mkdir -p results/part_d
make

T=8
C=64
N=4000
DAYS=200
REPEATS=5

PARTITIONS=(strips blocks interleaved)

RAW="results/part_d/d3_partitions_raw.csv"

echo "threads,chunks,partition,run,seconds" > "$RAW"

for PARTITION in "${PARTITIONS[@]}"; do

    for RUN in $(seq 1 "$REPEATS"); do

        echo "Testing $PARTITION, repetition $RUN"

        OUTPUT=$(./pandemic_sim \
            -N "$N" -D "$DAYS" -s 42 \
            -T "$T" -C "$C" \
            --pool --partition "$PARTITION" \
            --fixed-days --benchmark)

        ELAPSED_TIME=$(echo "$OUTPUT" |
            awk '/Simulation loop time:/ {print $4}')

        if [ -z "$ELAPSED_TIME" ]; then
            echo "Timing extraction failed"
            exit 1
        fi

        echo "$T,$C,$PARTITION,$RUN,$ELAPSED_TIME" >> "$RAW"

    done
done

echo "D3.2 completed"

