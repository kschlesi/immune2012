% 1D pathogen diffusion code w. all lymphocytes
% uses ode45 with @dyD1stoch

clear

global r_ sigma_ de_ f_ k_ days stepsize c gammas1D h_ ;

r_ = 3.3;
h_ = 10^-3;
sigma_ = 3;
de_ = 0.35;
k_ = 10^5;
f_ = 0.1;
c = 0.5;
b = 50;
N0density = 3;
days = 20;
stepsize = 0.1;
%mrate = 0.7; % per cell per day


% dimensions of 1D shape space
Pdim1 = 1000;
Ldim1 = 1000;

% gammas
gammas1D = zeros(Pdim1,Ldim1);
for i=1:Pdim1;
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end


% center and max amount of initial gaussian inoculation in shape space
x0 = 500;
Pmax0 = 10;
Pdiff0 = 12;

% setting initial conditions for P, N, E, M;
G0 = Pmax0.*ones(Pdiff0/2,1);
P0 = padarray(G0,Pdim1/2-Pdiff0/4,'both');

N0 = N0density.*ones(Ldim1,1);
E0 = zeros(Ldim1,1);
M0 = zeros(Ldim1,1);

% P0 = zeros(Pdim1,1); % initial gaussian distribution of pathogen
% for i=1:Pdim1;
%     P0(i) = Pmax0*exp(-1*((i-x0)^2)/(2*Pdiff0^2));
% end

% creating initial conditions vector

y0 = [P0;N0;E0;M0];

% integrating all diffeq in time
options = odeset('AbsTol',1e-3);
[ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),(0:stepsize:days),y0,options);
n_ts = size(ts_vec,1);

% create plotting functions

P_out = y_out(:,1:Pdim1);
N_out = y_out(:,Pdim1+1:Pdim1+Ldim1);
E_out = y_out(:,Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M_out = y_out(:,Pdim1+2*Ldim1+1:end);

% save results!!!
dlmwrite('Pdif11.txt',P_out);
dlmwrite('Ndif11.txt',N_out);
dlmwrite('Edif11.txt',E_out);
dlmwrite('Mdif11.txt',M_out);


% plot initial & final distributions
%     figure
%     plot((1:1:Pdim1),P0)
%     
%     figure
%     Pfin = squeeze(P_out(n_ts,:));
%     plot((1:1:Pdim1),Pfin)
    
    figure
    hold on
    surf(P_out,'MeshStyle','row')
    hold off
    axis([0 Pdim1 0 n_ts 0 max(P_out(:,x0))])
%    set(gca,'ZScale','log')

    figure
    hold on
    surf(N_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 N0density])

    figure
    hold on
    surf(E_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 5000])
    
    figure
    hold on
    surf(M_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 400])



    Ptot = sum(P_out,2);
    Ntot = sum(N_out,2);
    Etot = sum(E_out,2);
    Mtot = sum(M_out,2);
    figure
    semilogy(ts_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
    axis([0 days 1 10^8])
