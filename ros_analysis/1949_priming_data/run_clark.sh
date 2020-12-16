data1="prime_03_01-07-20-191923.txt"
data2="primer_12-20-19-150539_pta5nm_pta100nm_1949at5nm_1949_10nm_1949_100nm_1949_1uM_fls2_1uM_1949_1uM.txt"
data3="primer_1949_01-07-20-144240_pta5nm_1949at5nm_50nm1949_500nm1949_5um1949_19495umonly_fls2+5um1949+pta_1um1949.txt"
data4="primer_1949_changepta_01-07-20-201502_1um1949inall_pta5nm_pta+1949_pta10nm_pta+1949_pta100nm_pta+1949.txt"
data5="primer_pre_12-20-19-141256.txt"
data6="primer_pre_2_01-07-20-134942_1949.txt"


perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data1 -n 191923 -s 191923_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data2 -n 150539 -s 150539_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data3 -n 144240 -s 144240_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data4 -n 201502 -s 201502_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data5 -n 141256 -s 141256_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/1949_priming_data/ -r $data6 -n 134942 -s 134942_names.txt -tr -pn 1