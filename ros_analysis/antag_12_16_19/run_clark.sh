data1="antag_pta5uM_12-10-19-134034_pta_pa20_1410_1186_2009_2004_1949_5014.txt"
data2="antag_pta_14min_12-11-19-161257.txt"
data3="antag_5nM_5uM_12-20-19-132351_pta_pa20_1410_1186_2009_2004_1949_5014.txt"
data4="pre_antag_12-20-19-123102.txt"
data5="pre_antag_02-27-20-153934.txt"
data6="pta_antag_5nm_5uM_02-27-20-163305_pta_20_1410_5000_2009_2004_1949_5014.txt"

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data1 -n 134034 -s 134034_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data2 -n 161257 -s 161257_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data3 -n 132351 -s 132351_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data4 -n 123102 -s 123102_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data5 -n 153934 -s 153934_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/antag_12_16_19/ -r $data6 -n 163305 -s 163305_names.txt -tr -pn 1