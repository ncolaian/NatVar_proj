data1="col_fls2_bak1-5_10-30-19-160835_flg22_1391.txt"
data2="col_fls2_bak1_10-30-19-172709_flg22_1471_weird.txt"
data3="col_fls2_bak1_10-31-19-134623_flg22_1391.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak1_5_data_11_19/ -r $data1 -n 160835 -s 160835_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak1_5_data_11_19/ -r $data2 -n 172709 -s 172709_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/bak1_5_data_11_19/ -r $data3 -n 134623 -s 134623_names.txt -tr -pn 1