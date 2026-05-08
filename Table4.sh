#!/bin/bash

# 5h on a laptop (1 thread), mostly spent on the last three
PARAMS=$(cat <<-EOM
-n 128  -l 2  -L 2  --seed=0  --only-estimate
-n 128  -l 2  -L 1  --seed=0  --only-estimate
-n 128  -l 3  -L 1  --seed=0  --only-estimate
-n 128  -l 4  -L 1  --seed=0
-n 128  -l 4  -L 1  --seed=0  --only-estimate
-n 128  -l 5  -L 1  --seed=0  --only-estimate
-n 128  -l 6  -L 1  --seed=0  --only-estimate
-n 128  -l 7  -L 1  --seed=0  --only-estimate
-n 128  -l 8  -L 1  --seed=0  --only-estimate
-n 128  -l 9  -L 1  --seed=0  --only-estimate
-h
-n 192  -l 2  -L 2  --seed=0  --only-estimate
-n 192  -l 2  -L 1  --seed=0  --only-estimate
-n 192  -l 3  -L 1  --seed=0  --only-estimate
-n 192  -l 4  -L 1  --seed=0  --only-estimate
-n 192  -l 5  -L 1  --seed=0  --only-estimate
-n 192  -l 6  -L 1  --seed=0  --only-estimate
-n 192  -l 7  -L 1  --seed=0  --only-estimate
-n 192  -l 8  -L 1  --seed=0  --only-estimate
-n 192  -l 9  -L 1  --seed=0  --only-estimate
-h
-n 256  -l 3  -L 2  --seed=0  --only-estimate
-n 256  -l 3  -L 1  --seed=0  --only-estimate
-n 256  -l 4  -L 1  --seed=0  --only-estimate
-n 256  -l 5  -L 1  --seed=0  --only-estimate
-n 256  -l 6  -L 1  --seed=0  --only-estimate
-n 256  -l 7  -L 1  --seed=0  --only-estimate
-n 256  -l 8  -L 1  --seed=0  --only-estimate
-n 256  -l 9  -L 1  --seed=0  --only-estimate
-h
-n 128  -l 2  -L 127  --seed=0
-n 192  -l 2  -L 191  --seed=0
-n 256  -l 3  -L 382  --seed=0
-h
-n 128  -l 2  -L 15  --seed=0
-n 192  -l 2  -L 31  --seed=0
-n 256  -l 3  -L 94  --seed=0
EOM
)

IFS=$'\n'
for cmd in $PARAMS; do
	fname=$(echo $cmd | sed 's/  /_/g' | sed 's/[- ]//g')
	date
	IFS=$' '
	echo "CMD |$cmd| filename |$fname|"
	/usr/bin/time -v sage attack-aim2.py --no-progress-bars $cmd |& tee "logs/$fname.log"
	sleep 0.5
	IFS=$'\n'
done