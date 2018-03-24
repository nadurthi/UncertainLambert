clc
clear
filename_meta='../../cpp/examples/data/test_1.csv_meta';
filename_data='../../cpp/examples/data/test_1.csv_data';
[efmoutput,L,efmsize]=parse_cpp_struct(filename_meta, filename_data);


efmoutput=parse_cpp_efm(efmsize,efmoutput,L);


