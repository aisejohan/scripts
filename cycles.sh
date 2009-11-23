#! /bin/bash

change_data()
{
sed -i \
        -e "s@^#define d1\t.*@#define d1\t$d1@" \
        -e "s@^#define d2\t.*@#define d2\t$d2@" \
        -e "s@^#define d3\t.*@#define d3\t$d3@" \
        -e "s@^#define d4\t.*@#define d4\t$d4@" \
        -e "s@^#define d\t.*@#define d\t$d@" cycles/data.h
}

prepare_list_grobner()
{
	cat tmp/grobner-$d1-$d2-$d3-$d4-$d \
		| grep -v Coefficient \
		| grep -v Allocate \
		> tmp/list_grobner
}

make_cycles()
{
	cd cycles
	make input_pol
	cd ..
}

run_cycles()
{
	cat tmp/list_grobner | while read lijn
	do
		naam=$(echo "$lijn" | sed "s/ /-/g")
		grep -q "\[x - 2" tmp/$d1-$d2-$d3-$d4-$d/$naam
		if [ $? == 0 ]; then
			echo "This one: $naam"
			cd cycles
			echo $lijn | ./tester >> ../tmp/$d1-$d2-$d3-$d4-$d/$naam
			cd ..
		else
			echo "Not: $naam"
		fi
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
	change_data
	prepare_list_grobner
	make_cycles
	run_cycles
done
