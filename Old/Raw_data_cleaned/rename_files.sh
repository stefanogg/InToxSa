#!/bin/bash

for var in "$@"
do
	plate=$(echo $var | cut -d"/" -f2)
	exp_date=$(echo $var | cut -d"/" -f3)
	exp_date=${exp_date/_/""}
	exp_type="$(echo $var | cut -d"/" -f4)"
	new_name=${exp_type}_${exp_date}_${plate}
	filepath=$(dirname ${var})
  
	mkdir -p $exp_type
	mkdir -p ${exp_type}/xlsx
	mkdir -p ${exp_type}/csv
	
	echo "Copying "$var" to "${exp_type}/xlsx/${new_name}.xlsx
	cp $var ${exp_type}/xlsx/${new_name}.xlsx 

	echo "Converting to csv"
	echo $exp_type
	csv_file=${exp_type}/csv/${new_name}
	xlsx2csv ${exp_type}/xlsx/${new_name}.xlsx > ${csv_file}.csv 

	echo ""
	if [ $exp_type == "Growth_curves" ]
	then
	   echo "Removing headers on "$exp_type
	   drop_line=$(grep -n ",Time," Growth_curves/csv/${new_name}.csv |cut -d":" -f1)
	   echo "TEST"
	   echo $drop_line
	   tail -n +"$drop_line" ${csv_file}.csv > ${csv_file}_clean.csv
	fi
	   
	if [ $exp_type == "PI_kinetics" ]
	then
	   echo "Removing headers on "$exp_type
	   drop_line=$(grep -n "Well," PI_kinetics/csv/${new_name}.csv | cut -d":" -f1)
	   echo "TEST"
	   echo $drop_line
	   tail -n +"$drop_line" ${csv_file}.csv > ${csv_file}_clean.csv
	fi
done

