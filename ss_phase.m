% script for making various phase plots

clear

global r_ h_ sigma_ c beta_ chi_ Qstep x0 dh_ K_ ;
global b eps_ mu_ k_ Pdim1 Ldim1 delta_ Gamma_ muton ;

runnum = 2;
basecode = 'pldyn';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/']; %KONG datapath
% datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\' basecode '\']; %laptop datapath
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];
runnum = 1;
t2filename = [datapath 't' basecode num2str(runnum) '.txt'];
P2filename = [datapath 'P' basecode num2str(runnum) '.txt'];
L2filename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set global parameters, read in days from saved parameter file
params = setparams(bfilename);
days = params{end,2};    % total days run & saved in file

% data and time vector (read in from saved datafiles)
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Lplot = csvread(Lfilename);
t2plot = csvread(t2filename);
P2plot = csvread(P2filename);
L2plot = csvread(L2filename);

n_ts = size(tplot,1); % setting number of time steps
n_ts2 = size(t2plot,1);

% setting sub-mu_ populations to zero for plotting
for i=1:size(Lplot,1)
    for j=1:size(Lplot,2)
        if Lplot(i,j) < mu_
            Lplot(i,j) = 0;
        end
    end
end
for i=1:size(L2plot,1)
    for j=1:size(L2plot,2)
        if L2plot(i,j) < mu_
            L2plot(i,j) = 0;
        end
    end
end

% calculating total population vectors over time
    Ptot = sum(Pplot,2);
    Ltot = sum(Lplot,2);
    P2tot = sum(P2plot,2);
    L2tot = sum(L2plot,2);
    
% phase plot of total pathogen v. total lymphocytes
    figure
    loglog(P2tot,L2tot)
    title(['Phase-plot, ' basecode num2str(runnum)])
    xlabel('total pathogen')
    ylabel('total lymphocytes')
