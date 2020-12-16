data1="Data_08-29-18-115425_Pta_DA_146_264_166_1294_1810_1471_ROSAt.txt"
data2="data_08-29-18-125340_Pta_DA_1213_1544_1814_1184_Ps_Pa20.txt"
data3="data_08-29-18-160725_Pta_DA_255_133_102_99_321_118.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_29_18/ -r $data1 -n 115425 -s 115425_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_29_18/ -r $data2 -n 125340 -s 125340_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_29_18/ -r $data3 -n 160725 -s 160725_names.txt -tr -pn 1