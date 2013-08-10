clear

runnum = 1;
basecode = 'pldyn';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set parameters
params = setparams(bfilename);
for i=1:size(params,1)
    eval([char(params{i,1}) ' = ' num2str(params{i,2})]);
    eval([char(params{i,1}) 'units = char(params{i,3})']);
end
clear params;

