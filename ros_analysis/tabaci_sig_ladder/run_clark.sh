data1="Data_09-08-18-145522_flg22Pta_1_3_11_33_81_100_300_900_nM.txt"
data2="data_ladder_09-11-18-122236_1_3_10_20_33_81_100_300.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/tabaci_sig_ladder/ -r $data1 -n 145522 -s 145522_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/tabaci_sig_ladder/ -r $data2 -n 122236 -s 122236_names.txt -tr -pn 1