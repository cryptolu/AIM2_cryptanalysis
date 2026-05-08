#!/bin/bash

# n=128 15 IVs :       2^38 field ops

# n=256 70 IVs :  6GB  2^39 field ops

# n=128 11 IVs : 92GB  2^43 field ops
# n=192 20 IVs : 73GB  2^42 field ops
# n=256 52 IVs : 41GB  2^42 field ops

PARAMS=$(cat <<-EOM
-n 128  -l 2  -L 15  --seed=0
-n 192  -l 2  -L 24  --seed=0
-n 256  -l 3  -L 70  --seed=0
-h
-n 128  -l 2  -L 11  --seed=0
-n 192  -l 2  -L 20  --seed=0
-n 256  -l 3  -L 52  --seed=0
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