
data_pairs = [1 2 ; 1 6; 2 3; 2 4; 2 5; 6 7; 6 8];
my_gene = ["NR4A1" "NR4A2" "NR4A3" "ITBG1" "EHMT2" "LEF1" "AXIN2" "VIM" ...
    "ZEB1"  "CDH1"  "CDH3"];

root_path = pwd;
ndata = length(data_pairs);
for i = 1:ndata
    fprintf("Data %d %d \n", data_pairs(i,:));
    ipath = strcat( string(data_pairs(i,1)),"_");
    ipath = strcat( ipath, string(data_pairs(i,2)));
    cd(ipath);
    sce_file = strcat("sce", string(data_pairs(i,1)), string(data_pairs(i,2)));
    sce_file = strcat(sce_file, ".mat");
    load(sce_file);
    gene_boxplot(sce, my_gene);
    %fprintf("File %s \n", sce_file)
    cd(root_path);
end
