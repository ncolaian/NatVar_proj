data1="data_100nm_inhib_09-12-18-181238_pta_da_1297_59_20_pta1297_pta59_pta20.txt"
data2="data_100nm_inhib_09-12-18-191006_pta_da_1297_59_20_pta1297_pta59_pta20.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/inhib/ -r $data1 -n 181238 -s 181238_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/ros_analysis/inhib/ -r $data2 -n 191006 -s 191006_names.txt -tr -pn 1
