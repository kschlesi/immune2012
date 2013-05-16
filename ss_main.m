% 1D fitness/mutation infection integrator
% with cutoff implemented by looping calls to ode45

clear

global r_ h_ sigma_ de_ f_ k_ c b eps_ mu_ R_ dh_ K_ chi_ Qstep capon hsaton ;
global lambdas1D gammas1D tgone ntgone Nstep nrandon Gamma_ mrates delta_ ;

days = 10;      % number of days to run simulation
stepsize = 0.1; % size of steps at which to save data

% information about where to save data:
% this script will create 4 files whose names are defined here
runnum = ;
basecode = 'pldyn';
% datapath = ['C:\Documents and Settings\kimberly\My Documents\' ...
%     'Google Drive\immunedata\PL\' basecode '\']; %MOTHRA datapath
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\' basecode '\']; %laptop datapath
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% ensuring no overwrite of existing files
if isequal(exist(tfilename,'file'),2)
    error('Existing file; cannot overwrite!');
end

%%%%%%%%%%%%%%%%%%%%% setting necessary parameters %%%%%%%%%%%%%%%%%%%%%%%%
r_ = 3.3;           % pathogen mutation rate
h_ = 10^-5;         % pathogen killing
sigma_ = 3;         % naive recruitment
k_ = 10^5;          % pathogen saturation
c = 2;              % controls change of mutation prob. with distance
chi_ = 10;          % strength of mutation probability
Gamma_ = 4;         % naive influx
delta_ = 0.35;      % constant naive death rate
Qstep = 0.1;        % time-step for regenerating mutation matrix
Nstep = 5;          % time-step for regeneration naive cell distribution
b = 10;             % width of Gaussian affinity curve
beta_ = NaN;        % width of Gaussian fitness landscape
eps_ = 4;           % controls fall-off of fitness landscape at edges
mu_ = 1;            % minimum cell-per-site density
dh_ = 5*10^-7;      % coefficient of overall lymphocyte constraint
K_ = 10^10;         % pathogen carrying capacity
capon = 1;          % switches on/off pathogen carrying capacity
hsaton = 1;         % switches on/off lymphocyte constraint
nrandon = 0;        % switches on/off shuffling of naive cell distribution

% dimensions of 1D shape space
Pdim1 = 400;        % sites in pathogen shape space
Ldim1 = 400;        % sites in lymphocyte shape space
x0 = 6;             % location of original pathogen inoculation

% affinity and fitness information
gammas1D = zeros(Pdim1,Ldim1);  % matrix of affinities
lambdas1D = zeros(Pdim1,1);     % vector of pathogen fitnesses      
for i=1:Pdim1;
    lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end


%%%%%%%%%%%%%%%%%%%% setting initial configurations %%%%%%%%%%%%%%%%%%%%%%%
P0 = zeros(Pdim1,1);    % initial pathogen inoculation  
P0(4:8) = 3;    
% % initial gaussian distribution of pathogen
% P0 = zeros(Pdim1,1);
% for i=1:Pdim1;
%     P0(i) = Pmax0*exp(-1*((i-x0)^2)/(2*Pdiff0^2));
% end
L0density = Gamma_/delta_;          % initial naive cell mean density
%N0 = N0density.*ones(Ldim1,1);
L0 = unifrndpop(Ldim1,L0density,mu_); % random distribution of naive cells
R_ = Ldim1*L0density;   % total lymphocyte threshold, above which constraint applies


%%%%%%%%%%%%% writing parameters and init conditions to file %%%%%%%%%%%%%%
% saving/writing params to parameter file
b0 = [r_;h_;sigma_;k_;c;b;beta_;eps_;mu_;dh_;K_;R_;capon;hsaton;...
    Pdim1;Ldim1;x0;chi_;Qstep;Gamma_;Nstep;nrandon;delta_];
writeparams(bfilename,b0); % creates paramfile for run; returns error if file already exists

% creating & saving initial conditions vector
t0 = 0;
y0 = [P0;L0];

dlmwrite(tfilename,t0);
dlmwrite(Pfilename,transpose(P0));
dlmwrite(Lfilename,transpose(L0));


%%%%%%%%%%%%%%%%%%%%%%%% integrating diffeqs %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,mu_)); % sets 'events' option
tspan = (t0:stepsize:days); % timespan for solver                % to stop integration
n_ts = 1;       % counts total solver timesteps                  % and enforce cutoff at mu_
nstops = 0;     % counts number of interruptions
contin = 1;     % while loop parameter
tgone = 0;      % keeps track of most recent mutation matrix generation time 
ntgone = 0;     % keeps track of most recent naive cell redistribution time
mrates = eye(Pdim1);    % initial mutation matrix: no mutation

while (contin)
    
    % integrate until 'stopper' event...(or total days reached)
    % ('stopper.m' triggers an event whenever a population falls below mu_)
    [ts_vec,y_out,etimes,ytimes,indices] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),tspan,y0,options);
    
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
        cell2csv(afilename,tend,1); % appends cell line 'tend' to paramsfile
    end

end

%%%%%%%%%%%%%%%% plotting initial & final distributions %%%%%%%%%%%%%%%%%%%
    figure
    plot((1:1:Pdim1),P0)

    figure
    plot((1:1:Pdim1),P_out(end,:))