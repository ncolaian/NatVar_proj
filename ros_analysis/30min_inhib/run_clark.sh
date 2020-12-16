data1="Data_08-27-18-141714_Plate2_Pta_DA_457_94_236_468_289_470_ROS_At.txt"
data2="Data_08-27-18-171520_59_161_108_295_359_17_Da_Pta_ROSAt.txt"
data3="data_08-27-18-181200_Pta_DA_20_459_146_264_flg22Pa_flg20_ROSAt.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_27_18/ -r $data1 -n 141714 -s 141714_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_27_18/ -r $data2 -n 171520 -s 171520_names.txt -tr -pn 8

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_27_18/ -r $data3 -n 181200 -s 181200_names.txt -tr -pn 1