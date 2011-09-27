#! /bin/bash

cat tmp/list_done | while read line
do
	arr=($line)
	d1=${arr[0]}
	d2=${arr[1]}
	d3=${arr[2]}
	d4=${arr[3]}
	d=${arr[4]}
	grep "The characteristic polynomial is: " \
		tmp/$d1-$d2-$d3-$d4-$d/* | \
		sed 's@.*characteristic polynomial is: @@'
done

