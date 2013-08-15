% ss_seed
% given two previous runs PR1 and PR2, as well as a time t1, and 'days,' 
% this script will create a new run with the initial condition of PR1 at
% t1, using the paramfile of PR2, and run it until final time of 'days'+t1

clear

global mrates tgone ;

%%%%%%%%%%%% input information about seedfiles and newfile %%%%%%%%%%%%%%%%
PR1 = 'ftry1.2';  % run from which initial condition is drawn
PR2 = 'ftry1.2';  % run whose paramfile to use
t1 = 0;      % time in PR1 to use for initial condition; number or 'end'
days = 100;       % new days to append to file
stepsize = 0.1;  % size of steps at which to save

% new run files to be created
runnum = 1.4;
basecode = 'ftry';
isnew = 1;
datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
datapath1 = ['/Users/kimberly/Google Drive/immunedata/PL13/' deblank(char(PR1.*isletter(PR1))) '/'];
datapath2 = ['/Users/kimberly/Google Drive/immunedata/PL13/' deblank(char(PR2.*isletter(PR2))) '/'];

bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

b0filename = [datapath2 'b' PR2 '.txt'];
t0filename = [datapath1 't' PR1 '.txt'];
P0filename = [datapath1 'P' PR1 '.txt'];
L0filename = [datapath1 'L' PR1 '.txt'];
force_cont = 0;

% ensuring no overwrite of existing files
if isnew
    if isequal(exist(tfilename,'file'),2)
        error('Existing file; cannot overwrite!');
    end
end

% ensuring file existence
if (isequal(exist(t0filename,'file'),0)||isequal(exist(b0filename,'file'),0))
    error('At least one specified seedfile does not exist!');
end
if ~isnew
    if (isequal(exist(tfilename,'file'),0)||isequal(exist(bfilename,'file'),0))
        error('The file(s) to be continued do not exist!');
    end
end

% ensuring continuation reads from end and keeps old paramfile
if ~isnew
    if(~strcmp(tfilename,t0filename))
        error('File continuation appears to have been set in error!')
    end
    if(~strcmp(t1,'end') && ~strcmp(force_cont,'force_cont'))
        error(['File continuation must proceed from "end" -- to override,' ...
               'set force_cont = "force_cont"'])
    end
    if(~strcmp(bfilename,b0filename) && ~strcmp(force_cont,'force_cont'))
        error(['File continuation must use same paramfile -- to override,'...
               'set force_cont = "force_cont"'])
    end
end 


%%%%%%%%%%%%%% read in parameters and initial conditions %%%%%%%%%%%%%%%%%%

% set parameters, read in olddays
params = setparams(b0filename);
for i=1:size(params,1)
    if ~strcmp(char(params{i,1}),'days')
        eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
        eval([char(params{i,1}) 'units = char(params{i,3});']);
    end
end
olddays = t1;
if (strcmp(t1,'end'))
    olddays = params{end,2};    % days already run & saved in file
end
clear params;

% gammas & lambdas & mrates
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = Lambdas(eps_,Pdim1);     % vector of pathogen fitnesses      
for i=1:Pdim1;
    gammas1D(i,:) = Gammas([i,1],ones(Ldim1,1),1,b);
end
mrates = eye(Pdim1,Pdim1);
if (muton)
    mrates = Qmatrix(Pdim1,chi_);
end

% read in initial conditions (t1 from P/L0filename)
oldtimes = csvread(t0filename);
if(~strcmp(t1,'end'))
    t1in = find((oldtimes>t1));
    t0index = t1in(1)-1;
else
    t0index = size(oldtimes,1);
end
t0 = oldtimes(t0index);
P0 = transpose(csvread(P0filename,t0index-1,0,[t0index-1,0,t0index-1,Pdim1-1]));
L0 = transpose(csvread(L0filename,t0index-1,0,[t0index-1,0,t0index-1,Ldim1-1]));

% % modifying initial conditions vector (new infection?)
% P0_add = zeros(size(P0));
% P0_add(300:301) = 2;
% P0 = P0 + P0_add;


figure    % plot of P0 and L0 distributions at t0
hold on
hold all
plot((1:1:400),P0)
plot((1:1:400),L0)
title(['P0 and L0 seeding distributions at t = ' num2str(olddays) ' days'])
xlabel('location in shape space (site)')
ylabel('population (cells/\mul)')
legend('Pathogen','Lymphocytes')

figure
plot((1:1:400),lambdas1D)

%%%%%%%%%%%%% writing parameters and init conditions to file %%%%%%%%%%%%%%
% saving/writing params to parameter file
b0 = [r_;h_;sigma_;k_;b;eps_;mu_;dh_;K_;R_;capon;hsaton;...
    Pdim1;Ldim1;x0;chi_;Qstep;Gamma_;nrandon;delta_;muton;pinit];
if isnew
    writeparams(bfilename,b0); % creates paramfile for run; returns error if file already exists
end

% creating & saving initial conditions vector
y0 = [P0;L0];

if isnew
    dlmwrite(tfilename,t0);
    dlmwrite(Pfilename,transpose(P0));
    dlmwrite(Lfilename,transpose(L0));
end


%%%%%%%%%%%%%%%%%%%%%%%% integrating diffeqs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = (t0:stepsize:days+olddays);
options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,mu_));
n_ts = size(oldtimes,1);
contin = 1;
nstops = 0;
tgone = t0;
while (contin)

    % integrate until 'stopper' event...(or total days reached)
    % ('stopper.m' triggers an event whenever a population falls below mu_)
    [ts_vec,y_out,etimes,ytimes,indices] = ode45(@(t,y)ss_dy(t,y,b0,gammas1D,lambdas1D),...
        tspan,y0,options);

    % once integration is stopped...
    % add new internal steps to overall n_ts
    n_ts = n_ts + size(ts_vec,1)-1;
    nstops = nstops+1;

    % implement cutoff for site that fell below mu_
    if(indices)  % vector 'indices' contains site(s) that caused event
        for i=1:size(indices,1)
            if(y_out(end,indices(i,1))<=mu_)
                y_out(end,indices(i,1))=0;
            end
        end
    end
        
   % save & append y-output (leaving out old init condition)
    t_out = ts_vec(2:end);
    P_out = y_out(2:end,1:Pdim1);
    L_out = y_out(2:end,Pdim1+1:Pdim1+Ldim1);

    dlmwrite(tfilename,t_out,'-append');
    dlmwrite(Pfilename,P_out,'-append');
    dlmwrite(Lfilename,L_out,'-append');
    
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
       cell2csv(bfilename,tend,1); % appends cell line 'tend' to paramsfile
    end

end

%plot final distributions

figure    % plot of P0 and L0 distributions at days+olddays
hold on
hold all
plot((1:1:400),squeeze(P_out(end,:)))
plot((1:1:400),squeeze(L_out(end,:)))
title(['P0 and L0 distributions at t = ' num2str(ts_vec(end)) ' days'])
xlabel('location in shape space (site)')
ylabel('population (cells/\mul)')
legend('Pathogen','Lymphocytes')
