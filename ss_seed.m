% ss_seed
% given two previous runs PR1 and PR2, as well as a time t1, and 'days,' 
% this script will create a new run with the initial condition of PR1 at
% t1, using the paramfile of PR2, and run it until final time of 'days'+t1

clear
global r_ h_ sigma_ k_ c dh_ capon hsaton nrandon mrates Gamma_ delta_ muton ;
global b eps_ mu_ Pdim1 Ldim1 x0 chi_ Qstep Nstep tgone ntgone gammas1D lambdas1D;

%%%%%%%%%%%% input information about seedfiles and newfile %%%%%%%%%%%%%%%%
PR1 = 'pldyn1';  % run from which initial condition is drawn
PR2 = 'pldyn1';  % run whose paramfile to use
t1 = 0;          % time in PR1 to use for initial condition
days = 50;       % new days to append to file
stepsize = 0.1;  % size of steps at which to save

% new run files to be created
runnum = 4;
basecode = 'pldyn';
% datapath = ['C:\Documents and Settings\kimberly\My Documents\' ...
%     'Google Drive\immunedata\PL\']; %MOTHRA datapath
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\PL\']; %laptop datapath
bfilename = [datapath basecode '\b' basecode num2str(runnum) '.txt'];
tfilename = [datapath basecode '\t' basecode num2str(runnum) '.txt'];
Pfilename = [datapath basecode '\P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath basecode '\L' basecode num2str(runnum) '.txt'];

b0filename = [datapath deblank(char(PR2.*isletter(PR2))) '\b' PR2 '.txt'];
t0filename = [datapath deblank(char(PR1.*isletter(PR1))) '\t' PR1 '.txt'];
P0filename = [datapath deblank(char(PR1.*isletter(PR1))) '\P' PR1 '.txt'];
L0filename = [datapath deblank(char(PR1.*isletter(PR1))) '\L' PR1 '.txt'];

% ensuring no overwrite of existing files
if isequal(exist(tfilename,'file'),2)
    error('Existing file; cannot overwrite!');
end

% ensuring file existence
if (isequal(exist(t0filename,'file'),0)||isequal(exist(b0filename,'file'),0))
    error('At least one specified seedfile does not exist!');
end

% set parameters, read in olddays
params = setparams(b0filename);
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

mrates = eye(Pdim1,Pdim1);
if (muton)
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
end

%%%%%%%%%%%%% writing parameters and init conditions to file %%%%%%%%%%%%%%
% saving/writing params to parameter file
b0 = [r_;h_;sigma_;k_;c;b;beta_;eps_;mu_;dh_;K_;R_;capon;hsaton;...
    Pdim1;Ldim1;x0;chi_;Qstep;Gamma_;Nstep;nrandon;delta_;muton];
writeparams(bfilename,b0); % creates paramfile for run; returns error if file already exists

