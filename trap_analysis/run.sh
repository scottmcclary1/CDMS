#!/bin/bash 

bold=$(tput bold)
normal=$(tput sgr0)
base_dir=/N/dc2/scratch/scmcclar/school/graduate/spring2017/i590/i590_project
data=instrument_data_small
show_webpage=0

echo 
echo =============================================
echo ${bold}General Information${normal} 
echo =============================================
echo USER=$USER
date
echo HOSTNAME=$HOSTNAME

echo 
echo =============================================
echo ${bold}Move to Correct Location${normal} 
echo =============================================
cd $base_dir/trap_analysis
pwd 

echo 
echo =============================================
echo ${bold}Configure Environment${normal} 
echo =============================================
module load intel
ulimit -S -s 120000

echo 
echo =============================================
echo ${bold}Set OpenMP Environment Variables${normal}
echo =============================================
export OMP_STACKSIZE=100M
export OMP_NUM_THREADS=16
echo OMP_STACKSIZE=$OMP_STACKSIZE
echo OMP_NUM_THREADS=$OMP_STACKSIZE

echo
echo =============================================
echo ${bold}Build Executable${normal}
echo =============================================
make clean
make trap_omp
make clean

echo 
echo =============================================
echo ${bold}Move to Data Directory${normal}
echo =============================================
cd $base_dir/trap_analysis/$data
pwd

echo 
echo =============================================
echo ${bold}Run Application${normal}
echo =============================================
$base_dir/trap_analysis/./trap_omp < $base_dir/trap_analysis/input.txt
return_val=$?

echo 
echo =============================================
echo ${bold}Move Output${normal}
echo =============================================
mv -v data_Sep15_2014_* event_types_Sep15_2014_run08.txt peaks_Sep15_2014_run08.txt user_input_Sep15_2014.txt $base_dir/trap_analysis/output/

echo 
echo =============================================
echo ${bold}Move to Output Directory${normal}
echo =============================================
cd $base_dir/trap_analysis/output
pwd


echo 
echo =============================================
echo ${bold}Create Database Table\(s\)${normal}
echo =============================================
mysql -e "use scmcclar; drop table tbl_proj_data_filtered;"
mysql -e "CREATE TABLE scmcclar.tbl_proj_data_filtered ( "File" VARCHAR(12) NULL DEFAULT NULL , "Section_num" INT(1) NULL DEFAULT NULL , "End_Time_s" DOUBLE NULL DEFAULT NULL , "num_Points_Avgd" INT NULL DEFAULT NULL , "Freq_Avg" DOUBLE NULL DEFAULT NULL , "Freq_Slope" DOUBLE NULL DEFAULT NULL , "Freq_Rel_Slope" DOUBLE NULL DEFAULT NULL , "Freq_Sum_Sq" DOUBLE NULL DEFAULT NULL , "mz" DOUBLE NULL DEFAULT NULL , "mz_Std_Dev" DOUBLE NULL DEFAULT NULL , "Charge" DOUBLE NULL DEFAULT NULL , "Charge_Std_Dev" DOUBLE NULL DEFAULT NULL , "Mass" DOUBLE NULL DEFAULT NULL , "Harm2_Rel_Mag_Avg" DOUBLE NULL DEFAULT NULL , "Harm2_Rel_Mag_Std_Dev" DOUBLE NULL DEFAULT NULL);";

mysql -e "use scmcclar; drop table tbl_proj_event_types_single;"
mysql -e "CREATE TABLE scmcclar.tbl_proj_event_types_single ( "Single_Ion_File" VARCHAR(12) NULL DEFAULT NULL , "Section_num" INT(1) NULL DEFAULT NULL);";

mysql -e "use scmcclar; drop table tbl_proj_event_types_multiple;"
mysql -e "CREATE TABLE scmcclar.tbl_proj_event_types_multiple ( "Multiple_Ion_File" VARCHAR(12) NULL DEFAULT NULL , "Section_num" INT(1) NULL DEFAULT NULL);";

mysql -e "use scmcclar; drop table tbl_proj_event_types_no;"
mysql -e "CREATE TABLE scmcclar.tbl_proj_event_types_no ( "No_Ion_File" VARCHAR(12) NULL DEFAULT NULL , "Section_num" INT(1) NULL DEFAULT NULL);";

echo 
echo =============================================
echo ${bold}Ingest Output into Databases${normal}
echo =============================================
sed 's/ \{1,\}/,/g' data_Sep15_2014_filtered_run08.txt | tail -n +2 | sed '/^$/d' > ingest.txt
mysql -e "use scmcclar; LOAD DATA LOCAL INFILE 'ingest.txt' INTO TABLE tbl_proj_data_filtered FIELDS TERMINATED BY ',';";

sed 's/ \{1,\}/,/g' event_types_Sep15_2014_run08.txt | cut -d ',' -f 1-2 | tail -n +2 | sed '/^$/d' > ingest.txt
mysql -e "use scmcclar; LOAD DATA LOCAL INFILE 'ingest.txt' INTO TABLE tbl_proj_event_types_single FIELDS TERMINATED BY ',';";

sed 's/ \{1,\}/,/g' event_types_Sep15_2014_run08.txt | cut -d ',' -f 3-4 | tail -n +2 | sed '/^$/d' > ingest.txt
mysql -e "use scmcclar; LOAD DATA LOCAL INFILE 'ingest.txt' INTO TABLE tbl_proj_event_types_multiple FIELDS TERMINATED BY ',';";

sed 's/ \{1,\}/,/g' event_types_Sep15_2014_run08.txt | cut -d ',' -f 5-6 | tail -n +2 | sed '/^$/d' > ingest.txt
mysql -e "use scmcclar; LOAD DATA LOCAL INFILE 'ingest.txt' INTO TABLE tbl_proj_event_types_no FIELDS TERMINATED BY ',';";

rm ingest.txt

echo 
echo =============================================
echo ${bold}Show Small Part of Table\(s\)${normal}
echo =============================================
mysql -e "select File, Mass from scmcclar.tbl_proj_data_filtered limit 5;";
mysql -e "select * from scmcclar.tbl_proj_event_types_single limit 5;";
mysql -e "select * from scmcclar.tbl_proj_event_types_multiple limit 5;";
mysql -e "select * from scmcclar.tbl_proj_event_types_no limit 5;";

#select Single_Ion_File as File, tbl_proj_event_types_single.Section_num as Single_Ion_Section_num, tbl_proj_event_types_multiple.Section_num as Multiple_Ion_Section_num from tbl_proj_event_types_single join tbl_proj_event_types_multiple on tbl_proj_event_types_single.Single_Ion_File=tbl_proj_event_types_multiple.Multiple_Ion_File;


mysql -e "use scmcclar; select s.Single_Ion_File as File, s.Section_num as Single_Ion_Section_num, m.Section_num as Multiple_Ion_Section_num, n.Section_num as No_Ion_Section_num from tbl_proj_event_types_single as s left join tbl_proj_event_types_multiple as m on s.Single_Ion_File = m.Multiple_Ion_File left join tbl_proj_event_types_no as n on s.Single_Ion_File = n.No_Ion_File limit 10;";


#select s.Single_Ion_File as File, s.Section_num as Single_Ion_Section_num, m.Section_num as Multiple_Ion_Section_num, n.Section_num as No_Ion_Section_num from tbl_proj_event_types_single as s, tbl_proj_event_types_multiple m, tbl_proj_event_types_no as n where s.Single_Ion_File = m.Multiple_Ion_File and s.Single_Ion_File = n.No_Ion_File;

#select s.Single_Ion_File as File, s.Section_num as Single_Ion_Section_num, m.Section_num as Multiple_Ion_Section_num from tbl_proj_event_types_single as s left join tbl_proj_event_types_multiple as m on s.Single_Ion_File = m.Multiple_Ion_File limit 10;


if [ $show_webpage -eq "0" ]; then
    echo 
    echo =============================================
    echo ${bold}Open Webpage${normal}
    echo =============================================
    /N/u/scmcclar/Karst/local/src/nwjs-v0.15.1-linux-x64/nw $base_dir/webpage > /dev/null 2>&1
fi

echo 
echo =============================================
echo ${bold}Done${normal}
echo =============================================

exit $?
