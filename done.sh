#! /bin/bash

cat tmp/list | while read line
do
	arr=($line)
	d1=${arr[0]}
	d2=${arr[1]}
	d3=${arr[2]}
	d4=${arr[3]}
	d=${arr[4]}
	grep -q "$d1 $d2 $d3 $d4 $d" tmp/list_done
	if [ $? == 0 ]; then
		echo "This case: $d1 $d2 $d3 $d4 $d already done"
		exit 1
	fi
	echo $line >> tmp/list_done
done
