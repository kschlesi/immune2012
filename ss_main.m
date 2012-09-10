% simulation: pathogen diffused in 2D shape space
% interacts with lattice of lymphocytes

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
days = 20;
stepsize = 0.1;

% dimensions of 2D shape space
Ldim1 = 60;
Ldim2 = 48;
Pdim1 = 60; 
Pdim2 = 48;

% center and max amount of initial gaussian inoculation in shape space
x0 = [30,24];
Pmax0 = 4;
Pdiff0 = 5;

% setting initial conditions for P, N, E, M
%G0 = Gammas([20,20],zeros(40,40),Pmax0,Pdiff0);  % initial restricted gaussian distribution of pathogen
%P0 = padarray(G0,[Pdim1/2-20 Pdim2/2-20],'both');
P0 = Gammas(x0,zeros(Pdim1,Pdim2),Pmax0,Pdiff0);%+ones(Pdim1,Pdim2); % initial gaussian distribution of pathogen everywhere
N0 = N0density.*(ones(Ldim1,Ldim2)); % initial uniform distribution of naive cells
E0 = zeros(Ldim1,Ldim2);
M0 = zeros(Ldim1,Ldim2);

% option: restrict initial pathogen
% for i=1:Pdim1
%     for j=1:Pdim2
%         if (P0(i,j)<1)
%             P0(i,j)=0;
%         end
%     end
% end

% creating gamma library
gamma_lib = zeros(Pdim1,Pdim2,Ldim1,Ldim2);
for i = 1:Pdim1
    for j = 1:Pdim2
        gamma_lib(i,j,:,:) = Gammas([i,j],N0,1,b_);
    end
end

% %plotting test gamma
%         hold on
%         figure
%         surf(squeeze(gamma_lib(45,10,:,:)))
%         hold off

% creating vector of initial conditions
P0lin = reshape(P0,Pdim1*Pdim2,1);
N0lin = reshape(N0,Ldim1*Ldim2,1);
E0lin = reshape(E0,Ldim1*Ldim2,1);
M0lin = reshape(M0,Ldim1*Ldim2,1);
y0 = [P0lin;N0lin;E0lin;M0lin];

% plotting initial conditions
    hold on
    figure
    surf(P0)
    hold off

 % integrating diffeqs
 [ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Pdim2,Ldim1,Ldim2),(0:stepsize:days),y0);
 n_ts = size(ts_vec,1);
 
% setting final values
P_outlin = y_out(:,1:Pdim1*Pdim2); % 2D matrices of cells per site (linear) per ts
N_outlin = y_out(:,Pdim1*Pdim2+1:Pdim1*Pdim2+Ldim1*Ldim2);
E_outlin = y_out(:,Pdim1*Pdim2+Ldim1*Ldim2+1:Pdim1*Pdim2+2*Ldim1*Ldim2);
M_outlin = y_out(:,Pdim1*Pdim2+2*Ldim1*Ldim2+1:end);

dlmwrite('Pout17.txt',P_outlin);
dlmwrite('Nout17.txt',N_outlin);
dlmwrite('Eout17.txt',E_outlin);
dlmwrite('Mout17.txt',M_outlin);

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
        
