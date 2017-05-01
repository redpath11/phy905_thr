#!/bin/bash

TEMP=550
while [ $TEMP -lt 601 ]; do
	echo "Starting Temperature: ${TEMP}K"
	#Save statistics.txt for each simulation
	FILEMV="Stats_${TEMP}.txt"
	echo "Saved statistics file is ${FILEMV}"

	./moldyn.exe 5 $TEMP 5.26
	mv statistics.txt ${FILEMV}

	let TEMP=TEMP+5
done
