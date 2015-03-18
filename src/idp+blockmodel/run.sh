#!/bin/sh

for i in `seq 4 10`; do
#	echo $i
	java Main friend $i 0 > log/$i.log &
done

java Main friend 20 0 > log/20.log &
java Main friend 30 0 > log/30.log &
java Main friend 40 0 > log/40.log &
java Main friend 50 0 > log/50.log &

