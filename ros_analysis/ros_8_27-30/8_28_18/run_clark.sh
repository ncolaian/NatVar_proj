data1="Data_08-28-18-132713_Pta_DA_1814_1184_255_133_102_99_ROSAt.txt"
data2="Data_08-28-18-151250_Pta_DA_321_118_457_94_236_468_ROSAt.txt"
data3="Data_08-28-18-160734_Pta_DA_289_470_17_359_295_108_ROSAt.txt"
data4="data_08-28-18_pta_da_161_59_20_459_pa_pa20.txt"
data5="data_08-28-18_pta_da_166_1294_1810_1471_1213_1544.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_28_18/ -r $data1 -n 132713 -s 132713_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_28_18/ -r $data2 -n 151250 -s 151250_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_28_18/ -r $data3 -n 160734 -s 160734_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_28_18/ -r $data4 -n 161_59 -s 161_59_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/8_28_18/ -r $data5 -n 166_1294 -s 166_1294_names.txt -tr -pn 1