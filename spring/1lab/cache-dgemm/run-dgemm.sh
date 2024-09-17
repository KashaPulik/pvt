#!/bin/sh

# taskset --cpu-list 0 ./dgemm
perf stat -e cache-misses taskset --cpu-list 0 ./dgemm
