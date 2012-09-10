% simulation: pathogen diffused in shape space
% interacts with lattice of 

clear

global b_ r_ h_ sigma_ k_ de_ f_ gamma_lib days stepsize;

N0density = 3; 
r_ = 3.3;
h_ = 10^-5;
sigma_ = 3;
k_ = 10^5;
de_ = 0.35;
f_ = 0.1;
b_ = 5;
days = 10;
stepsize = 0.1;
olddays = 2;
oldss = 0.1;
oldt_vec = (0:oldss:olddays);
old_ts = size(oldt_vec,2);

% dimensions of 2D shape space
Ldim1 = 60;
Ldim2 = 48;
Pdim1 = 60; 
Pdim2 = 48;

% %plotting test gamma
%         hold on
%         figure
%         surf(squeeze(gamma_lib(45,10,:,:)))
%         hold off

% creating vector of initial conditions
P0in = csvread('Pout17.txt');
N0in = csvread('Nout17.txt');
E0in = csvread('Eout17.txt');
M0in = csvread('Mout17.txt');
P0lin = shiftdim(P0in(old_ts,:),1);
N0lin = shiftdim(N0in(old_ts,:),1);
E0lin = shiftdim(E0in(old_ts,:),1);
M0lin = shiftdim(M0in(old_ts,:),1);

y0 = [P0lin;N0lin;E0lin;M0lin];


% plotting initial conditions
    P0 = reshape(P0lin,Pdim1,Pdim2);
    hold on
    figure
    surf(P0)
    hold off
    
    N0 = reshape(N0lin,Ldim1,Ldim2);
    hold on
    figure
    surf(N0)
    hold off
   
 % creating gamma library
gamma_lib = zeros(Pdim1,Pdim2,Ldim1,Ldim2);
for i = 1:Pdim1
    for j = 1:Pdim2
        gamma_lib(i,j,:,:) = Gammas([i,j],N0,1,b_);
    end
end

 % integrating diffeqs
 options = odeset('AbsTol',1e-3);
 [ts_vec,y_out] = ode45(@(t,y)sssDiffusedy_sf1(t,y,Pdim1,Pdim2,Ldim1,Ldim2),(0:stepsize:days),y0,options);
 n_ts = size(ts_vec,1);
 
% setting final values
P_outlin = y_out(:,1:Pdim1*Pdim2); % 2D matrices of cells per site (linear) per ts
N_outlin = y_out(:,Pdim1*Pdim2+1:Pdim1*Pdim2+Ldim1*Ldim2);
E_outlin = y_out(:,Pdim1*Pdim2+Ldim1*Ldim2+1:Pdim1*Pdim2+2*Ldim1*Ldim2);
M_outlin = y_out(:,Pdim1*Pdim2+2*Ldim1*Ldim2+1:end);

dlmwrite('Pnewbie.txt',P_outlin);
dlmwrite('Nnewbie.txt',N_outlin);
dlmwrite('Enewbie.txt',E_outlin);
dlmwrite('Mnewbie.txt',M_outlin);

concat('Pout17.txt','Pnewbie.txt');
concat('Nout17.txt','Nnewbie.txt');
concat('Eout17.txt','Enewbie.txt');
concat('Mout17.txt','Mnewbie.txt');

P_out = reshape(P_outlin,n_ts,Pdim1,Pdim2); % 3D matrices of cells per row per col per ts
N_out = reshape(N_outlin,n_ts,Ldim1,Ldim2);
E_out = reshape(E_outlin,n_ts,Ldim1,Ldim2);
M_out = reshape(M_outlin,n_ts,Ldim1,Ldim2);

P_tot = sum(P_outlin,2); % 1D vectors of total cells per ts
N_tot = sum(N_outlin,2); 
E_tot = sum(E_outlin,2);
M_tot = sum(M_outlin,2);

% plotting immune system response to pathogen
    figure
        semilogy(ts_vec,P_tot,ts_vec,(M_tot + N_tot + E_tot))%,ts_vec,M_tot,ts_vec,N_tot,ts_vec,E_tot);
        axis([0 days 1 10^10])

% plotting final P, N, M, and E populations
        hold on
        figure
        surf(squeeze(P_out(n_ts,:,:)))
        hold off

        hold on
        figure
        surf(squeeze(N_out(n_ts,:,:)))
        hold off

        hold on
        figure
        surf(squeeze(E_out(n_ts,:,:)))
        hold off

        hold on
        figure
        surf(squeeze(M_out(n_ts,:,:)))
        hold off
        
