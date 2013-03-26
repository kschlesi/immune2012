function fixparams()

% copypaste in before y0 created
% add NaN params that are not defined in code
% figure out totaldays and set t0 as such

afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
a0 = [r_;h_;sigma_;de_;k_;f_;c;b;beta_;eps_;mu_;dh_;K_;R_;capon;hsaton;Pdim1;Ldim1;x0];
writeparams(afilename,a0) % creates paramfile for run; returns error if file already exists
t0 = ;
tend = cell(1,3);
tend{1,1} = 'days';
tend{1,2} = t0;
tend{1,3} = 'days';
cell2csv(afilename,tend,1); % appends cell line 'tend' to paramsfile