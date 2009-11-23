#! /bin/bash

prepare_list_grobner()
{
	cat tmp/grobner-$d1-$d2-$d3-$d4-$d \
		| grep -v Coefficient \
		| grep -v Allocate \
		> tmp/list_grobner
}

run_gp()
{
	cat tmp/list_grobner | while read lijn
	do
		naam=$(echo "$lijn" | sed "s/ /-/g")
		gp -q tmp/$d1-$d2-$d3-$d4-$d/$naam do_it.gp >> tmp/$d1-$d2-$d3-$d4-$d/$naam
	done
}

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
	prepare_list_grobner
	run_gp
done
