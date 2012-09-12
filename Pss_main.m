% one-D diffusion equation test (discrete)

clear

global r_ days stepsize mrate D gammas ;

r_ = 3.3;
D = 30;
days = 5;
stepsize = 0.1;
mrate = 0.7; % per cell per day

% dimensions of 1D shape space
Pdim1 = 200;

% center and max amount of initial gaussian inoculation in shape space
x0 = 100;
Pmax0 = 10;
Pdiff0 = 12;

% setting initial conditions for P
G0 = Pmax0.*ones(Pdiff0/2,1);
P0 = padarray(G0,Pdim1/2-Pdiff0/4,'both');

% P0 = zeros(Pdim1,1); % initial gaussian distribution of pathogen
% for i=1:Pdim1;
%     P0(i) = Pmax0*exp(-1*((i-x0)^2)/(2*Pdiff0^2));
% end

% creating initial conditions vector


% integrating all diffeq in time
options = odeset('AbsTol',1e-3);
[ts_vec,P_out] = ode45(@(t,y)Pss_dy(t,y,Pdim1),[0:stepsize:days],P0,options);
n_ts = size(ts_vec,1);

% save results!!!
dlmwrite('P1diftest.txt',P_out);

% plot initial & final P distributions
    figure
    plot([1:1:Pdim1],P0)
    
    figure
    Pfin = squeeze(P_out(n_ts,:));
    plot([1:1:Pdim1],Pfin)
    
    figure
    hold on
    surf(P_out,'MeshStyle','row')
    hold off
    axis([0 Pdim1 0 n_ts 0 max(P_out(:,x0))])
%    set(gca,'ZScale','log')
