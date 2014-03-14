% 1D fitness/mutation infection integrator
% with cutoff implemented by looping calls to ode45

clear

global mrates ;

days = 100;       % number of days to run simulation
stepsize = 0.1; % size of steps at which to save data

% information about where to save data:
% this script will create 4 files whose names are defined here
runnum = 4.7;
basecode = 'plos';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% ensuring no overwrite of existing files
if isequal(exist(tfilename,'file'),2)
    error('Existing file; cannot overwrite!');
end

%%%%%%%%%%%%%%%%%%%%% setting necessary parameters %%%%%%%%%%%%%%%%%%%%%%%%
r_ = 3;             % pathogen mutation rate
h_ = 10^-5;         % pathogen killing
sigma_ = 3;         % naive recruitment
k_ = 10^5;          % pathogen saturation
chi_ = 22;          % strength of mutation probability (chi_=0: no mutation)
Gamma_ = 1;         % naive influx
delta_ = 0.33;      % constant naive death rate
pinit = 10;         % initial dose of pathogen
b = 23.9;             % width of Gaussian affinity curve
eps_ = 0;           % controls fall-off of fitness landscape at edges
mu_ = 1;            % minimum cell-per-site density
dh_ = 5e-7;         % coefficient of overall lymphocyte constraint
K_ = 10^10;         % pathogen carrying capacity
capon = 1;          % switches on/off pathogen carrying capacity
hsaton = 1;         % switches on/off lymphocyte constraint
nrandon = 0;        % switches on/off stochastic naive distribution
spliton = 0;        % switched on/off mutation term split from growth term

% dimensions of 1D shape space
Pdim1 = 400;        % sites in pathogen shape space
Ldim1 = 400;        % sites in lymphocyte shape space
x0 = 6;             % location of original pathogen inoculation

% affinity and fitness information
gammas1D = zeros(Pdim1,Ldim1);   % matrix of affinities
lambdas1D = Lambdas(eps_,Pdim1); % vector of pathogen fitnesses      
for i=1:Pdim1;
    gammas1D(i,:) = Gammas([i,1],ones(Ldim1,1),1,b);
end
mrates = Qmatrix(Pdim1,chi_,spliton);    % initial mutation matrix


%%%%%%%%%%%%%%%%%%%% setting initial configurations %%%%%%%%%%%%%%%%%%%%%%%
P0 = zeros(Pdim1,1);    % initial pathogen inoculation  
P0(x0) = pinit;    

Qprime = 1;  % ratio of R to L_tot* (i.e. R = L* times n times Qprime)
L0density = Gamma_/(delta_ - dh_*(1-Qprime));          % initial naive cell mean density
if (nrandon)
    L0 = unifrndpop(Ldim1,L0density,mu_); % random distribution of naive cells
else
    L0 = L0density.*ones(Ldim1,1)-1;      % uniform distribution of naive cells
end
R_ = Ldim1*L0density*Qprime;   % total lymphocyte threshold, above which constraint applies


%%%%%%%%%%%%% writing parameters and init conditions to file %%%%%%%%%%%%%%
% saving/writing params to parameter file
b0 = [r_;h_;sigma_;k_;b;eps_;mu_;dh_;K_;R_;capon;hsaton;...
    Pdim1;Ldim1;x0;chi_;Gamma_;nrandon;delta_;spliton;pinit];
writeparams(bfilename,b0); % creates paramfile for run; returns error if file already exists

% creating & saving initial conditions vector
t0 = 0;
y0 = [P0;L0];

dlmwrite(tfilename,t0);
dlmwrite(Pfilename,P0');
dlmwrite(Lfilename,L0');


%%%%%%%%%%%%%%%%%%%%%%%% integrating diffeqs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,mu_)); % sets 'events' option
tspan = (t0:stepsize:days); % timespan for solver                % to stop integration
n_ts = 1;       % counts total solver timesteps                  % and enforce cutoff at mu_
nstops = 0;     % counts number of interruptions
contin = 1;     % while loop parameter
while (contin)
    
    % integrate until 'stopper' event...(or total days reached)
    % ('stopper.m' triggers an event whenever a population falls below mu_)
    [ts_vec,y_out,etimes,ytimes,indices] = ode45(@(t,y)ss_dy(t,y,b0,gammas1D,lambdas1D),...
        tspan,y0,options);
    
    % once integration is stopped...
    % add new internal steps to overall n_ts
    n_ts = n_ts + size(ts_vec,1)-1;
    nstops = nstops + 1;
    
    % implement cutoff for site that fell below mu_
    if(indices)  % vector 'indices' contains site(s) that caused event
        for i=1:size(indices,1)
            if(y_out(end,indices(i,1))<=mu_)
                y_out(end,indices(i,1))=0;
            end
        end
    end
        
    % save & append ode45 output (leaving out old init condition)
    t_out = ts_vec(2:end);
    P_out = y_out(2:end,1:Pdim1);        
    L_out = y_out(2:end,Pdim1+1:Pdim1+Ldim1);

    dlmwrite(tfilename,t_out,'-append');
    dlmwrite(Pfilename,P_out,'-append');
    dlmwrite(Lfilename,L_out,'-append');
    
    % set new initial conditions for resuming while loop integration
    t0 = ts_vec(end);
    tspan = (t0:stepsize:days);
    y0 = y_out(end,:);
    
    if(t0>=days-stepsize)  % if within one stepsize of the end
        tspan = [t0,days]; % reset tspan to default steps
    end
    
    if(t0>=days) % check stopping condition; if end reached:
        contin = 0;                         % end while loop
        tend = cell(1,3);                   % create cell of final 'days'
        tend{1,1} = 'days';                 % and write it to paramfile
        tend{1,2} = t0;
        tend{1,3} = 'days';
        cell2csv(bfilename,tend,1); % appends cell line 'tend' to paramsfile
    end

end

%%%%%%%%%%%%%%%% plotting initial & final distributions %%%%%%%%%%%%%%%%%%%

figure    % plot of P0 and L0 distributions at days
plot((1:1:400),P_out(end,:)')
hold on
hold all
plot((1:1:400),L_out(end,:)')
title(['P0 and L0 distributions at t = ' num2str(ts_vec(end)) ' days'])
xlabel('location in shape space (site)')
ylabel('population (cells/\mul)')
legend('Pathogen','Lymphocytes')
