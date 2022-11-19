#!/bin/bash
./benchmark.sh lcsh-big -n 400 -r 20 -t"1 10 20 40 80" -p"-v"
mv lcsh-big-bench.log lcsh-big-bench_exact.log
./benchmark.sh lcsh-ramaeu -n 400 -r 20 -t"1 10 20 40 80" -p"-v"
mv lcsh-ramaeu-bench.log lcsh-ramaeu-bench_exact.log
./benchmark.sh homo-muso -n 400 -r 20 -t"1 10 20 40 80" -p"-v"
mv homo-muso-bench.log homo-muso-bench_exact.log
./benchmark.sh dmela-scere -n 400 -r 20 -t"1 10 20 40 80" -p"-v"
mv dmela-scere-bench.log dmela-scere-bench_exact.log

