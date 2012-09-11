% test for diffused antigen mutation in shape space

clear

global r_ days stepsize std Qmat odiag;

r_ = 3.3;
days = 5;
stepsize = 0.1;
std = sqrt(0.1);

% dimensions of 2D shape space
Pdim1 = 20; 
Pdim2 = 20;

% center and max amount of initial gaussian inoculation in shape space
x0 = [10,10];
Pmax0 = 10;
Pdiff0 = 2;

% Q-matrix must be Pdim1 x Pdim2 x Pdim1 x Pdim2
% preallocate:
Qmat = zeros(Pdim1,Pdim2,Pdim1,Pdim2);
odiag = zeros(4,1);

% setting initial conditions for P
G0 = Gammas([5,5],zeros(10,10),Pmax0,Pdiff0);  % initial restricted gaussian distribution of pathogen
P0 = padarray(G0,[Pdim1/2-5 Pdim2/2-5],'both');
%P0 = Gammas(x0,zeros(Pdim1,Pdim2),Pmax0,Pdiff0);  % initial gaussian distribution of pathogen
y0 = reshape(P0,Pdim1*Pdim2,1);

% integrating P diffeq in time
options = odeset('AbsTol',1e-3);
[ts_vec,y_out] = ode45(@(t,y)stochastic(t,y,Pdim1,Pdim2),[0:stepsize:days],y0,options);
n_ts = size(ts_vec,1);
P_out = reshape(y_out,n_ts,Pdim1,Pdim2);

% save results!!!
dlmwrite('Pmut5c.txt',y_out);

% plot initial & final P distributions
    hold on
    figure
    surf(P0)
    hold off

    figure
    hold on
    Pfin = squeeze(P_out(n_ts,:,:));
    surf(Pfin)
    hold off