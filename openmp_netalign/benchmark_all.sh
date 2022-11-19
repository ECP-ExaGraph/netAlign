#!/bin/bash

die () {
    echo >&2 "$@"
    exit 1
}

PROB_DIR="../data-mtx"
PROBLEM=$1
OUTPUT=$PROBLEM-bench.log
NITER=50

# clear the "$1" option 
shift $(( 1 ))

# default, don't build
BUILD=
BATCH=0
EPARAMS=

nthreads="1 2 4 8 16 32"

while getopts ":bn:qt:r:p:" opt; do
  case $opt in
    b)
      BUILD=1
      ;;
    n)
      NITER=$OPTARG
      ;;
    q)
      OUTPUT=/dev/null
      ;;
    t)
      nthreads=$OPTARG
      ;;
    r) 
      BATCH=$OPTARG
      ;;
    p) 
      EPARAMS=$OPTARG
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      exit 1
      ;;
  esac
done

OUTPUT=$PROBLEM-r$BATCH-bench.log

[ -e $PROB_DIR/$PROBLEM-A.mtx ] || die "Problem file $PROB_DIR/$PROBLEM-A.mtx does not exist"

if [ ! -z "$BUILD" ]; then
    echo "Building fresh copy ..."
    make clean
    make
fi

export OMP_NESTED=1
export KMP_AFFINITY=granularity=fine,scatter,0,1

PARAMS="-n $NITER -r $BATCH -f $EPARAMS"

echo | tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo "Benchmarking performance for problem"| tee --append $OUTPUT
echo "  " $PROBLEM| tee --append $OUTPUT
echo "on "`date`| tee --append $OUTPUT
echo "  with KMP_AFFINITY=$KMP_AFFINITY"| tee --append $OUTPUT
echo "  with numactl --interleave=all "| tee --append $OUTPUT
echo "  with PARAMS=$PARAMS"| tee --append $OUTPUT
echo "  with NTHREADS=$nthreads"| tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo | tee --append $OUTPUT

for nt in $nthreads; do
  echo 
  echo "Threads $nt" | tee --append $OUTPUT
  export OMP_NUM_THREADS=$nt
  if [ $nt ==  1 ]; then
		command="numactl --membind=all ./netAlign $PROB_DIR/$PROBLEM $PARAMS"
  else
		command="numactl --interleave=all ./netAlign $PROB_DIR/$PROBLEM $PARAMS"
  fi
  echo "Running $command" | tee --append $OUTPUT 
  
  # run command, putting everything into output, 
  # but also sending it to the scree (/dev/tty)
  # and then finding the "Solve Time" 
  
  STIME=`$command | tee --append $OUTPUT | tee /dev/tty | grep "Solve Time: " | cut -f3 -d' '`
  times[$nt]=$STIME
done

ONETIME=${times[1]}

echo | tee --append $OUTPUT
echo "Final Report"| tee --append $OUTPUT
for nt in $nthreads; do
  divexp="$ONETIME / ${times[$nt]}"
  speedup=$(awk "BEGIN{print $divexp}")
  printf "  Threads %4i : %7.2f s.  %4.1fx speedup\\n" $nt ${times[$nt]} $speedup| tee --append $OUTPUT
done

export KMP_AFFINITY=granularity=fine,compact,0,1

echo | tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo "Benchmarking performance for problem"| tee --append $OUTPUT
echo "  " $PROBLEM| tee --append $OUTPUT
echo "on "`date`| tee --append $OUTPUT
echo "  with KMP_AFFINITY=$KMP_AFFINITY"| tee --append $OUTPUT
echo "  with numactl --interleve=all"| tee --append $OUTPUT
echo "  with PARAMS=$PARAMS"| tee --append $OUTPUT
echo "  with NTHREADS=$nthreads"| tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo | tee --append $OUTPUT

for nt in $nthreads; do
  echo 
  if [ $nt != 1 ]; then
  	echo "Threads $nt" | tee --append $OUTPUT
  	export OMP_NUM_THREADS=$nt
	command="numactl --interleave=all ./netAlign $PROB_DIR/$PROBLEM $PARAMS"
  echo "Running $command" | tee --append $OUTPUT 
  
  # run command, putting everything into output, 
  # but also sending it to the scree (/dev/tty)
  # and then finding the "Solve Time" 
  
  STIME=`$command | tee --append $OUTPUT | tee /dev/tty | grep "Solve Time: " | cut -f3 -d' '`
  times[$nt]=$STIME
  fi
done


echo | tee --append $OUTPUT
echo "Final Report"| tee --append $OUTPUT
for nt in $nthreads; do
  divexp="$ONETIME / ${times[$nt]}"
  speedup=$(awk "BEGIN{print $divexp}")
  printf "  Threads %4i : %7.2f s.  %4.1fx speedup\\n" $nt ${times[$nt]} $speedup| tee --append $OUTPUT
done

export KMP_AFFINITY=granularity=fine,scatter,0,1

echo | tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo "Benchmarking performance for problem"| tee --append $OUTPUT
echo "  " $PROBLEM| tee --append $OUTPUT
echo "on "`date`| tee --append $OUTPUT
echo "  with KMP_AFFINITY=$KMP_AFFINITY"| tee --append $OUTPUT
echo "  with numactl --membind=all"| tee --append $OUTPUT
echo "  with PARAMS=$PARAMS"| tee --append $OUTPUT
echo "  with NTHREADS=$nthreads"| tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo | tee --append $OUTPUT

for nt in $nthreads; do
  echo 
  if [ $nt != 1 ]; then
  	echo "Threads $nt" | tee --append $OUTPUT
  	export OMP_NUM_THREADS=$nt
  	command="numactl --membind=all ./netAlign $PROB_DIR/$PROBLEM $PARAMS"
  	echo "Running $command" | tee --append $OUTPUT 
  
  # run command, putting everything into output, 
  # but also sending it to the scree (/dev/tty)
  # and then finding the "Solve Time" 
  
  	STIME=`$command | tee --append $OUTPUT | tee /dev/tty | grep "Solve Time: " | cut -f3 -d' '`
  	times[$nt]=$STIME
  fi
done


echo| tee --append $OUTPUT 
echo "Final Report"| tee --append $OUTPUT
for nt in $nthreads; do
  divexp="$ONETIME / ${times[$nt]}"
  speedup=$(awk "BEGIN{print $divexp}")
  printf "  Threads %4i : %7.2f s.  %4.1fx speedup\\n" $nt ${times[$nt]} $speedup| tee --append $OUTPUT
done

export KMP_AFFINITY=granularity=fine,compact,0,1

echo | tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo "Benchmarking performance for problem"| tee --append $OUTPUT
echo "  " $PROBLEM| tee --append $OUTPUT
echo "on "`date`| tee --append $OUTPUT
echo "  with KMP_AFFINITY=$KMP_AFFINITY"| tee --append $OUTPUT
echo "  with numactl --membind=all"| tee --append $OUTPUT
echo "  with PARAMS=$PARAMS"| tee --append $OUTPUT
echo "  with NTHREADS=$nthreads"| tee --append $OUTPUT
echo "-------------------------------------"| tee --append $OUTPUT
echo | tee --append $OUTPUT

for nt in $nthreads; do
  echo 
  if [ $nt != 1 ]; then
  	echo "Threads $nt" | tee --append $OUTPUT
  	export OMP_NUM_THREADS=$nt
  	command="numactl --membind=all ./netAlign $PROB_DIR/$PROBLEM $PARAMS"
  	echo "Running $command" | tee --append $OUTPUT 
  
  # run command, putting everything into output, 
  # but also sending it to the scree (/dev/tty)
  # and then finding the "Solve Time" 
  
  	STIME=`$command | tee --append $OUTPUT | tee /dev/tty | grep "Solve Time: " | cut -f3 -d' '`
  	times[$nt]=$STIME
  fi
done


echo | tee --append $OUTPUT
echo "Final Report"| tee --append $OUTPUT
for nt in $nthreads; do
  divexp="$ONETIME / ${times[$nt]}"
  speedup=$(awk "BEGIN{print $divexp}")
  printf "  Threads %4i : %7.2f s.  %4.1fx speedup\\n" $nt ${times[$nt]} $speedup| tee --append $OUTPUT
done
