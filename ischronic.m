%CHRONICITY

function [ischronic,strains] = ischronic(basecode,runnum)

datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
%Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set parameters
params = setparams(bfilename);
for i=1:size(params,1)
    eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
    eval([char(params{i,1}) 'units = char(params{i,3});']);
end
clear params;

% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
%Lplot = csvread(Lfilename);

changes = diff((Pplot>=mu_)')';
strains = max(sum((changes>0),2),sum((changes<0),2));

strains = strains + (Pplot(:,1)>=mu_).*(Pplot(:,end)>=mu_);

figure
plot(tplot,strains)

ischronic=0;

end