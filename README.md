# scGeatoobox_addson
Complementary functions for scGeatoobox from CaiLab
Pulling the code via MATLAB terminal

tic;
disp('Getting scGEAToolbox_addson...');
unzip('https://github.com/sromero1992/scGEAToolbox_addson/archive/main.zip');
addpath('./scGEAToolbox_addson-main');
toc;
savepath(fullfile(userpath,'pathdef.m'));
