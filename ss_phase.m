% script for making various phase plots

clear

global r_ h_ sigma_ c beta_ chi_ Qstep x0 dh_ K_ ;
global b eps_ mu_ k_ Pdim1 Ldim1 delta_ Gamma_ muton ;

pdays = 50;
runnum = 2;
basecode = 'pldyn';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/']; %KONG datapath
% datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\' basecode '\']; %laptop datapath
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];
runnum = 1;
b2filename = [datapath 'b' basecode num2str(runnum) '.txt'];
t2filename = [datapath 't' basecode num2str(runnum) '.txt'];
P2filename = [datapath 'P' basecode num2str(runnum) '.txt'];
L2filename = [datapath 'L' basecode num2str(runnum) '.txt'];
runnum = 5;
b3filename = [datapath 'b' basecode num2str(runnum) '.txt'];
t3filename = [datapath 't' basecode num2str(runnum) '.txt'];
P3filename = [datapath 'P' basecode num2str(runnum) '.txt'];
L3filename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set global parameters, read in days from saved parameter file
params = setparams(bfilename);
days = params{end,2}(1,1);    % total days run & saved in file
params = setparams(b2filename);
days2 = params{end,2}(1,1);
params = setparams(b3filename);
days3 = params{end,2}(1,1);

% data and time vector (read in from saved datafiles)
tplot = csvread(tfilename);
if(pdays < days)
    t1in = find((tplot > pdays));
    t0index = t1in(1)-1;
else
    t0index = size(tplot,1);
end
Pplot = (csvread(Pfilename,0,0,[0,0,t0index-1,Pdim1-1]));
Lplot = (csvread(Lfilename,0,0,[0,0,t0index-1,Ldim1-1]));
t2plot = csvread(t2filename);
if(pdays < days2)
    t1in = find((t2plot > pdays));
    t0index = t1in(1)-1;
else
    t0index = size(t2plot,1);
end
P2plot = (csvread(P2filename,0,0,[0,0,t0index-1,Pdim1-1]));
L2plot = (csvread(L2filename,0,0,[0,0,t0index-1,Ldim1-1]));
t3plot = csvread(t3filename);
if(pdays < days3)
    t1in = find((t3plot > pdays));
    t0index = t1in(1)-1;
else
    t0index = size(t3plot,1);
end
P3plot = (csvread(P3filename,0,0,[0,0,t0index-1,Pdim1-1]));
L3plot = (csvread(L3filename,0,0,[0,0,t0index-1,Ldim1-1]));

n_ts = size(tplot,1); % setting number of time steps
n_ts2 = size(t2plot,1);
n_ts3 = size(t3plot,1);

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
for i=1:size(L3plot,1)
    for j=1:size(L3plot,2)
        if L3plot(i,j) < mu_
            L3plot(i,j) = 0;
        end
    end
end


% calculating total population vectors over time
    Ptot = sum(Pplot,2);
    Ltot = sum(Lplot,2);
    P2tot = sum(P2plot,2);
    L2tot = sum(L2plot,2);
    P3tot = sum(P3plot,2);
    L3tot = sum(L3plot,2);
    
% phase plot of total pathogen v. total lymphocytes
    figure
    loglog(Ptot,Ltot,'--')
    hold on
    hold all
    loglog(P3tot,L3tot,'-.')
    loglog(P2tot,L2tot)
    title('Phase-plot dynamics')
    xlabel('total pathogen density')
    ylabel('total lymphocyte density')
    legend('infection (i)','infection (ii)','infecton(iii)')
