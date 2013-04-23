% 1D fitness/mutation infection integrator
% with cutoff implemented by looping calls to ode45

clear

global r_ h_ sigma_ de_ f_ k_ c b eps_ mu_ R_ dh_ K_ chi_ Qstep capon hsaton ;
global lambdas1D gammas1D tgone mrates ;

days = 10;
stepsize = 0.1; % size of steps at which to save

runnum = 2;
basecode = 'naive';
%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\' basecode '\'];%NEW laptop Gdrive
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\'; %laptop datapath
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% ensuring no overwrite
if isequal(exist(tfilename,'file'),2)
    error('Existing file; cannot overwrite!');
end

%%%%%%%%%%%%%%% setting necessary parameters %%%%%%%%%%%%%%%%%%%
r_ = 3.3;
h_ = 10^-5;
sigma_ = 3;
de_ = 0.35;
k_ = 10^5;
f_ = 0.1;
c = 2;
chi_ = 10;
Qstep = 0.1;
b = 10;
beta_ = NaN; 
eps_ = 4; 
mu_ = 1;
dh_ = 5*10^-7;
K_ = 10^10;
capon = 1;
hsaton = 1;

% dimensions of 1D shape space
Pdim1 = 400;
Ldim1 = 400;
x0 = 6;

% gammas & lambdas
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end


%%%%%%%%%%%%%%%% setting initial configurations %%%%%%%%%%%%%%%%%%%
P0 = zeros(Pdim1,1);
%P0(x0-2:x0+2) = 3;
P0(4:8) = 3;
% % initial gaussian distribution of pathogen
% P0 = zeros(Pdim1,1);
% for i=1:Pdim1;
%     P0(i) = Pmax0*exp(-1*((i-x0)^2)/(2*Pdiff0^2));
% end
N0density = 3;
%N0 = N0density.*ones(Ldim1,1);
N0 = unifrndpop(Ldim1,N0density,mu_);
E0 = zeros(Ldim1,1);
M0 = zeros(Ldim1,1);
R_ = Ldim1*N0density;


%%%%%%%%%%% writing parameters and init conditions to file %%%%%%%%%%%
% saving/writing params to paramfile
a0 = [r_;h_;sigma_;de_;k_;f_;c;b;beta_;eps_;mu_;dh_;K_;R_;capon;hsaton;Pdim1;Ldim1;x0;chi_;Qstep];
writeparams(afilename,a0); % creates paramfile for run; returns error if file already exists

% creating initial conditions vector
t0 = 0;
y0 = [P0;N0;E0;M0];

dlmwrite(tfilename,t0);
dlmwrite(Pfilename,transpose(P0));
dlmwrite(Nfilename,transpose(N0));
dlmwrite(Efilename,transpose(E0));
dlmwrite(Mfilename,transpose(M0));    


%%%%%%%%%%%%%%%%%%% integrating diffeqs %%%%%%%%%%%%%%%%%%%%%%%%%%%
options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,mu_));
tspan = (t0:stepsize:days);
n_ts = 1;
nstops = 0;
contin = 1;
tgone = 0;
mrates = eye(Pdim1);
while (contin)
    
    % integrate until 'stopper' event...(or days)
    [ts_vec,y_out,etimes,ytimes,indices] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),tspan,y0,options);
        
    % add new internal steps to overall n_ts
    size(ts_vec,1)
    n_ts = n_ts + size(ts_vec,1)-1;
    nstops = nstops + 1;
    
    % implement one-cell cutoff for PNEM that caused event ONLY(!)
    if(indices)
        for i=1:size(indices,1)
            if(y_out(end,indices(i,1))<=mu_)
                y_out(end,indices(i,1))=0;
            end
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
    
    % set new initial conditions for resuming while loop integration
    t0 = ts_vec(end);
    tspan = (t0:stepsize:days);
    y0 = y_out(end,:);
    
    if(t0>=days-stepsize)  % if within one stepsize of final days
        tspan = [t0,days]; % reset tspan to default steps
    end
    
    if(t0>=days) % check stopping condition; if days reached:
        contin = 0;                         % end while loop
        tend = cell(1,3);                   % create cell of final 'days'
        tend{1,1} = 'days';                 % and write it to paramfile
        tend{1,2} = t0;
        tend{1,3} = 'days';
        cell2csv(afilename,tend,1); % appends cell line 'tend' to paramsfile
    end

end

%%%%%%%%%%%%%% plotting initial & final distributions %%%%%%%%%%%%%%%%%%%
    figure
    plot((1:1:Pdim1),P0)
    
    figure
    plot((1:1:Pdim1),P_out(end,:))