export base_path=$PWD
echo $base_path;
cd LHScans; mkdir plots;
cd $base_path; 
mkdir combine_files;
mkdir datacard; cd datacard; mkdir datacard_2016; mkdir datacard_2017; mkdir datacard_2018;
cd $base_path;
mv prepareCards.py datacard;
mkdir plots;
cd templates; mkdir 2016; mkdir 2017; mkdir 2018;
cd $base_path;
mkdir scripts;
cd $base_path;
