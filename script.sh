#! /bin/bash

echo "7 11 16 25 91 10 5 3 6 3" > tmp/list

change_data()
{
sed -i \
        -e "s@^#define d1\t.*@#define d1\t$d1@" \
        -e "s@^#define d2\t.*@#define d2\t$d2@" \
        -e "s@^#define d3\t.*@#define d3\t$d3@" \
        -e "s@^#define d4\t.*@#define d4\t$d4@" \
        -e "s@^#define d\t.*@#define d\t$d@" new/data.h
sed -i \
        -e "s@^#define d1\t.*@#define d1\t$d1@" \
        -e "s@^#define d2\t.*@#define d2\t$d2@" \
        -e "s@^#define d3\t.*@#define d3\t$d3@" \
        -e "s@^#define d4\t.*@#define d4\t$d4@" \
        -e "s@^#define d\t.*@#define d\t$d@" grobner/data.h
}

make_grobner()
{
	cd grobner
	make make_list
	cd ..
}

run_grobner()
{
	cd grobner
	./tester | tee ../tmp/grobner-$d1-$d2-$d3-$d4-$d
	cd ..
}

prepare_list_grobner()
{
	cat tmp/grobner-$d1-$d2-$d3-$d4-$d \
		| grep -v Coefficient \
		| grep -v Allocate \
		> tmp/list_grobner
}

make_new()
{
	cd new
	make input_pol
	cd ..
}

run_new()
{
	mkdir tmp/$d1-$d2-$d3-$d4-$d
	cd new
	cat ../tmp/list_grobner | while read lijn
	do
		naam=$(echo "$lijn" | sed "s/ /-/g")
		echo $lijn | ./tester | tee ../tmp/$d1-$d2-$d3-$d4-$d/$naam
	done
	cd ..
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
	make_grobner
	run_grobner
	prepare_list_grobner
	make_new
	run_new
done

