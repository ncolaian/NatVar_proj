data1="fe_last_peps_10-15-19-174413_cpa_fepa_cpa20_fepa20_c5006_fe5006_c5013_fe5013.txt"
data2="fe_last_peps_10-15-19-174420_cpta_fepta_cda_feda_c161_fe161_c1256_fe1256.txt"
data3="fe_last_peps_10-17-19-160736_cpa_fepa_c5007_fe5007.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/fls2efr_10_16_19/ -r $data1 -n 174413 -s 174413_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/fls2efr_10_16_19/ -r $data2 -n 174420 -s 174420_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/fls2efr_10_16_19/ -r $data3 -n 160736 -s 160736_names.txt -tr -pn 1
