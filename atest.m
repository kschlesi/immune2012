disp('RUNNING!!');

disp('testing new Rfull threshold thingy...');
Rfull = 10^6;

runnum = 1.4;
basecode = 'plos';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set parameters
days = 0;
params = setparams(bfilename);
for i=1:size(params,1)
    eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
    eval([char(params{i,1}) 'units = char(params{i,3});']);
end
clear params;

% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Lplot = csvread(Lfilename);

n_ts = size(tplot,1);

% update paramfile with last 'days' value saved if ode45 was interrupted
if ~days
    tend = cell(1,3);                   
    tend{1,1} = 'days';                
    tend{1,2} = tplot(end);
    tend{1,3} = 'days';
    cell2csv(bfilename,tend,1); 
    clear tend;
    days = tplot(end);
end    
    
Pplot = Pplot.*(Pplot>=mu_);
Lplot = Lplot.*(Lplot>=mu_);
Ptot = sum(Pplot,2);
Ltot = sum(Lplot,2);

nhidden = Rfull-n;
R = n*Gamma_/delta_;
factor = delta_/(5-1);
decaytermN = factor*(1-(Ltot+nhidden*Lplot(:,n-2*x0))/Rfull);
decaytermL = dh_*R*(1-Ltot/R);
figure
plot(decaytermN,decaytermL)
xlabel('new decay');
ylabel('old decay');
figure
plot(tplot,decaytermN,tplot,decaytermL)
legend('new decay','old decay');