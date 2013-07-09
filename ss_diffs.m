clear

global r_ h_ sigma_ k_ c dh_ K_ R_ capon hsaton nrandon mrates Gamma_ delta_ muton beta_;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

frametimes=205;

runnum = 1;
basecode = 'pldyn';
Pdim1 = 400;
Ldim1 = 400;

ss_min = 1;
ss_max = 100;
t_min = 4;
t_max = 400;

datapath = ['/Users/kimberly/Google Drive/immunedata/PL/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

params = setparams(bfilename);
days = params{end,2};

timevec = csvread(tfilename);

%find tstep_max
tinrange = find(timevec>t_max);
if numel(tinrange)==0
    tstep_max = size(timevec,1);
else
    tstep_max = tinrange(1)-1;
end

%find tstep_min
tinrange = find(timevec<t_min);
if numel(tinrange)==0
    tstep_min = 1;
else
    tstep_min = tinrange(end)+1;
end
clear tinrange;

% bound ss_max and ss_min
if(ss_min<1)
    ss_min = 1;
end
if(ss_max>Pdim1||ss_max>Ldim1)
    ss_max = min(Pdim1,Ldim1);
end
%disp([t_min t_max timevec(tstep_min) timevec(tstep_max)])
%disp([ss_max-ss_min+1 tstep_max-tstep_min+1])


% read in desired rectangle
P = (csvread(Pfilename,tstep_min-1,ss_min-1,[tstep_min-1,ss_min-1,tstep_max-1,ss_max-1]));
%L = transpose(csvread(Lfilename,tstep_min-1,ss_min-1,[tstep_min-1,ss_min-1,tstep_max-1,ss_max-1]));
%disp([size(P) size(L)]);

Xaxis = timevec(tstep_min:tstep_max);
Yaxis = transpose(1:1:ss_max);

figure
contourf(Xaxis,Yaxis,transpose(P>0),'EdgeColor','none');

figure
contourf(Xaxis(2:end),Yaxis,transpose(diff(P>0)>0),'EdgeColor','none');
