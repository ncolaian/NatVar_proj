data1="data_inhib_5nm_30min_10-04-18-181839_pta+pta_dw+pta_59+pta_59+dw_1297+pta_1297+dw_20+pta_20+dw.txt"
data2="data_inhib_5nm_30min_t1_10-04-18-120356_pta+pta_dw+pta_da+pta_da+dw_1292+pta_1292+dw_20+pta_20+dw.txt"
data3="data_inhib_5nm_30min_t2_10-04-18-132556_pts+pta_dw+pta_da+pta_da+dw_1292+pta_1292+dw_20+pta_20+dw.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/30min_inhib/10_04_18/ -r $data1 -n 181839 -s 181839_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/30min_inhib/10_04_18/ -r $data2 -n 120356 -s 120356_names.txt -tr -pn 8

perl /Users/nicholascolaianni/Documents/dangl_lab/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/Ros_analysis/30min_inhib/10_04_18/ -r $data3 -n 132556 -s 132556_names.txt -tr -pn 1