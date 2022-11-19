#!/bin/bash
#make clean
#make
export KMP_STACKSIZE=32m
export OMP_NESTED=1
>homoMP
echo "Threads 1">>homoMP
export OMP_NUM_THREADS=1
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 2">>homoMP
export OMP_NUM_THREADS=2
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 4">>homoMP
export OMP_NUM_THREADS=4
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 8">>homoMP
export OMP_NUM_THREADS=8
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 16">>homoMP
export OMP_NUM_THREADS=16
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 32">>homoMP
export OMP_NUM_THREADS=32
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
echo "Threads 64">>homoMP
export OMP_NUM_THREADS=64
./netAlign A.mtx B.mtx L.mtx 1 2 0.999 1000 0 3 >>homoMP
>dmelaMP
echo "Threads 1">>dmelaMP
export OMP_NUM_THREADS=1
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 2">>dmelaMP
export OMP_NUM_THREADS=2
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 4">>dmelaMP
export OMP_NUM_THREADS=4
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 8">>dmelaMP
export OMP_NUM_THREADS=8
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 16">>dmelaMP
export OMP_NUM_THREADS=16
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 32">>dmelaMP
export OMP_NUM_THREADS=32
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
echo "Threads 64">>dmelaMP
export OMP_NUM_THREADS=64
./netAlign A1.mtx B1.mtx L1.mtx 1 2 0.999 1000 0 3 >>dmelaMP
