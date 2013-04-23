% ss_main_cont

clear

global r_ h_ sigma_ de_ f_ k_ c dh_ capon hsaton nrandon mrates Gamma_ ;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

days = 2;        % new days to append to file
stepsize = 0.1;  % size of steps at which to save

% file to which new days will be appended
runnum = 2;
basecode = 'naive';
%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\' basecode '\'];%NEW laptop Gdrive
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\'; %laptop
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% ensuring file existence
if isequal(exist(tfilename,'file'),0)
    error('File to be continued does not exist!');
end

% set parameters, read in olddays
params = setparams(afilename);
olddays = params{end,2};    % days already run & saved in file

% gammas & lambdas & mrates
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

mrates = zeros(Pdim1,Pdim1);
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/(i-j)^c)/chi_;
            mrates(j,i) = mrates(i,j);
        end
        iloss = sum(mrates(i,:));
        mrates(i,i) = 1-iloss;
    end
iloss = sum(mrates(1,:));
mrates(1,1) = 1-iloss;

% initial inoculation in shape space
t0in = csvread(tfilename);
P0in = csvread(Pfilename);
N0in = csvread(Nfilename);
E0in = csvread(Efilename);
M0in = csvread(Mfilename);
old_ts = size(t0in,1);
P0 = shiftdim(P0in(old_ts,:),1);
N0 = shiftdim(N0in(old_ts,:),1);
E0 = shiftdim(E0in(old_ts,:),1);
M0 = shiftdim(M0in(old_ts,:),1);


% creating initial conditions vector
t0 = t0in(end);
y0 = [P0;N0;E0;M0];

% integrating diffeqs in time with a FOR LOOP
tspan = (t0:stepsize:days+olddays);
options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,Pdim1));
n_ts = old_ts;
contin = 1;
nstops = 0;
tgone = t0;
ntgone = t0;
while (contin)
    
    % integrate until jth event... (or days)
    [ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),tspan,y0,options);

    % add new internal steps to overall n_ts
    size(ts_vec,1)
    n_ts = n_ts + size(ts_vec,1)-1;
    nstops = nstops+1;
%     n_ts
%     nstops

    % implement one-cell cutoff for all P
    for pcount=1:Pdim1
        if(y_out(end,pcount)<=mu_)
            y_out(end,pcount)=0;
        end
    end
        
   % save & append y-output (leaving out old init condition)
    t_out = ts_vec(2:end);
    P_out = y_out(2:end,1:Pdim1);
    N_out = y_out(2:end,Pdim1+1:Pdim1+Ldim1);
    E_out = y_out(2:end,Pdim1+Ldim1+1:Pdim1+2*Ldim1);
    M_out = y_out(2:end,Pdim1+2*Ldim1+1:end);

    dlmwrite(tfilename,t_out,'-append');
    dlmwrite(Pfilename,P_out,'-append');
    dlmwrite(Nfilename,N_out,'-append');
    dlmwrite(Efilename,E_out,'-append');
    dlmwrite(Mfilename,M_out,'-append'); 
    
    % set new initial conditions
    t0 = ts_vec(end);
    y0 = y_out(end,:);
    tspan = (t0:stepsize:days+olddays);
        
    if(t0>=days+olddays-stepsize)
        tspan = [t0,days+olddays];
    end
    
    if (t0>=days+olddays)
       contin = 0; 
       tend = cell(1,3);
       tend{1,1} = 'days';
       tend{1,2} = t0;
       tend{1,3} = 'days';
       cell2csv(afilename,tend,1); % appends cell line 'tend' to paramsfile
    end

end


% plot initial & final distributions
    figure
    plot((1:1:Pdim1),P0)
    
    figure
    Pfin = squeeze(P_out(end,:));
    plot((1:1:Pdim1),Pfin)
    
%     Ptot = sum(P_out,2);
%     Ntot = sum(N_out,2);
%     Etot = sum(E_out,2);
%     Mtot = sum(M_out,2);
%     figure
%     semilogy(ts_vec(2:end),Ptot,ts_vec(2:end),Ntot+Mtot+Etot)
