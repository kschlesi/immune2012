function tlengths()

clear

% global r_ h_ sigma_ de_ f_ c beta_ chi_ x0 ;
% global b eps_ mu_ k_ Pdim1 Ldim1 ;

runnum = 3;
basecode = 'edge';
%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';%laptop
%datapath = 'C:\Users\Kimberly\Desktop\immune2012_data\'; %M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
% Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
% Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
% Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
% Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% set parameters, read in days
params = setparams(afilename);
days = params{end,2};    % total days run & saved in file

tplot = csvread(tfilename);
% Pplot = csvread(Pfilename);
% Nplot = csvread(Nfilename);
% Eplot = csvread(Efilename);
% Mplot = csvread(Mfilename);

% create vector of timestep lengths in order
tsteplengths = zeros(size(tplot));
for i=2:size(tplot,1)
    tsteplengths(i) = tplot(i)-tplot(i-1);
end

mean(tsteplengths)
std(tsteplengths)
days/size(tplot,1)

figure
scatter((1:1:size(tplot,1)),tsteplengths)
figure
hist(tsteplengths)