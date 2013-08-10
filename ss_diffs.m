clear

frametimes=205;

runnum = 2;
basecode = 'pldyn';
Pdim1 = 0;
Ldim1 = 0;

ss_min = 1;
ss_max = 100;
t_min = 4;
t_max = 400;

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
