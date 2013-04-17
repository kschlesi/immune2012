% script for making various phase plots

clear

global r_ h_ sigma_ de_ f_ c beta_ chi_ Qstep x0 dh_ K_ ;
global b eps_ mu_ k_ Pdim1 Ldim1 ;

runnum = 6;
basecode = 'edge';
%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\' basecode '\'];%NEW laptop Gdrive
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';%OLD laptop
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %OLD M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% set global parameters, read in days from saved parameter file
params = setparams(afilename);
days = params{end,2};    % total days run & saved in file

% data and time vector (read in from saved datafiles)
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(tplot,1); % setting number of time steps

% setting sub-mu_ populations to zero for plotting
for i=1:size(Nplot,1)
    for j=1:size(Nplot,2)
        if Nplot(i,j) < mu_
            Nplot(i,j) = 0;
        end
        if Eplot(i,j) < mu_
            Eplot(i,j) = 0;
        end
        if Mplot(i,j) < mu_
            Mplot(i,j) = 0;
        end
    end
end

% calculating total population vectors over time
    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    Ltot = Ntot + Mtot + Etot;
    
% phase plot of total pathogen v. total lymphocytes
    figure
    loglog(Ptot,Ltot)
    title(['Phase-plot, ' basecode num2str(runnum)])
    xlabel('total pathogen')
    ylabel('total lymphocytes')
