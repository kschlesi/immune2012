% ss_main_cont

clear

%global r_ h_ sigma_ de_ f_ k_ c ;
global b beta_ mu_;

days = 250;      % new days to append to file
stepsize = 0.1;  % size of steps at which to save
olddays = 50;    % days already run & saved in file
oldss = 0.1;

% file to which new days will be appended
datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\';
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';
tfilename = [datapath 'tbound18.txt'];
Pfilename = [datapath 'Pbound18.txt'];
Nfilename = [datapath 'Nbound18.txt'];
Efilename = [datapath 'Ebound18.txt'];
Mfilename = [datapath 'Mbound18.txt'];

% ensuring file existence
if isequal(exist(tfilename),0)
    error('File to be continued does not exist!');
end

% dimensions of 1D shape space
Pdim1 = 400;
Ldim1 = 400;
x0 = 200;

% gammas & lambdas
p_ = (1-exp(-1*((Pdim1)^2)/(8*beta_^2)))^(-1);
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - p_*(1-exp(-1*((i-x0)^2)/(2*beta_^2)));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% initial inoculation in shape space
t0in = csvread(tfilename);
P0in = csvread(Pfilename);
N0in = csvread(Nfilename);
E0in = csvread(Efilename);
M0in = csvread(Mfilename);
old_ts = size(P0in,1);
P0 = shiftdim(P0in(old_ts,:),1);
N0 = shiftdim(N0in(old_ts,:),1);
E0 = shiftdim(E0in(old_ts,:),1);
M0 = shiftdim(M0in(old_ts,:),1);


% creating initial conditions vector
t0 = t0in(end);
y0 = [P0;N0;E0;M0];

% plotting initial conitions
    % figure
    % semilogy((1:Pdim1+3*Ldim1),(shiftdim(y0)))

% integrating diffeqs in time with a FOR LOOP
tspan = (t0:stepsize:days+olddays);
options = odeset('AbsTol',1e-3,'Events',@(t,y)stopper(t,y,Pdim1));
n_ts = old_ts;
contin = 1;
nstops = 0;
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
    end

end


% plot initial & final distributions
    figure
    plot((1:1:Pdim1),P0)
    
    figure
    Pfin = squeeze(P_out(end,:));
    plot((1:1:Pdim1),Pfin)
    
    Ptot = sum(P_out,2);
    Ntot = sum(N_out,2);
    Etot = sum(E_out,2);
    Mtot = sum(M_out,2);
    figure
    semilogy(ts_vec(2:end),Ptot,ts_vec(2:end),Ntot+Mtot+Etot)
