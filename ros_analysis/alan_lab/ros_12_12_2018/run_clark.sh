data1="Data_12-12-18-152245_pta_ptada_pa_pa20_1297_447_20_470.txt"
data2="Data_12-12-18-162647_pta_1nm_3nm_7nm_10nm_33nm_81nm_100nm_300nm.txt"
data3="Data_12-13-18-155932_pta_ptada_pa_pa20_1297_1447_20_470.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/alan_lab/ros_12_12_2018/ -r $data1 -n 152245 -s 152245_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/alan_lab/ros_12_12_2018/ -r $data2 -n 162647 -s 162647_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/alan_lab/ros_12_12_2018/ -r $data3 -n 155932 -s 155932_names.txt -tr -pn 1