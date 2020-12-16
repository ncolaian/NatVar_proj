data1="otu5_antag_09-09-20-184838_m_pa205uM_pv2V5uM_potAntag5uM_ptEvade5uM_m_m_1949100nM.txt"
data2="otu5_pre_09-10-20-190612_m_pa205uM_pvAnt_potAntag_potEvader_m_m.txt"
data3="otu5_pre_09-10-20-190613_pta5nMx5_potAntag_potEvader_100nM.txt"
data4="otu_post_09-10-20-124933_pta5nMx5_Pant100nM_Pevade100nM_pta5nM.txt"
data5="post_otu5_09-19-20-175715_ptax5_pe_a.txt"
data6="post_otu5_09-19-20-192310_ptax5_pe_a.txt"
data7="pre_otu5_09-19-20-164930_m_20_e_a_aV_m_m.txt"
data8="pre_otu5_09-19-20-181725_m_20_e_a_aV_m_m.txt"


perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data1 -n 184838 -s 184838_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data2 -n 190612 -s 190612_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data3 -n 190613 -s 190613_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data4 -n 124933 -s 124933_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data5 -n 175715 -s 175715_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data6 -n 192310 -s 192310_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data7 -n 164930 -s 164930_names.txt -tr -pn 1

perl /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/clark_matlab_code/clark_runner.pl -o /Users/nicholascolaianni/Documents/dangl_lab/nat_variants_proj/ros_analysis/otu_5_ros/ -r $data8 -n 181725 -s 181725_names.txt -tr -pn 1