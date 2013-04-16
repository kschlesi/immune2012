clear

global r_ h_ sigma_ de_ f_ c beta_ chi_ Qstep x0 capon hsaton;
global b eps_ mu_ k_ Pdim1 Ldim1 dh_ K_ ;

runnum = 16;
basecode = 'qstep';
datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';%laptop
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% set parameters, read in days
params = setparams(afilename);
days = params{end,2};    % total days run & saved in file

% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(tplot,1);

for i=1:size(Nplot,1)
    for j=1:size(Nplot,2)
        if Nplot(i,j) < 1
            Nplot(i,j) = 0;
        end
        if Eplot(i,j) < 1
            Eplot(i,j) = 0;
        end
        if Mplot(i,j) < 1
            Mplot(i,j) = 0;
        end
    end
end

figure
plot((1:1:n_ts),Eplot(:,50))