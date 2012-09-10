% shape space simulation 1
% Sean's original; calls ShSpSimdy
% 2D, no x-shapespace

clear

global b_ r_ h_ sigma_ k_ de_ f_ gammas_ ;

N0density = 3; 
P0 = 10; 
r_ = 3.3;
h_ = 10^-5;
sigma_ = 2;
k_ = 10^5;
de_ = 0.35;
f_ = 0.1;
b_ = 5;
days = 30;

% dimensions of 2D shape space
dim1 = 60;
dim2 = 50;

% initial location of inoculation in shape space
x1 = [30,25];

% setting initial conditions for N, E, M
N0 = N0density.*(ones(dim1,dim2))-(N0density-0.01).*padarray(ones(dim1-40,dim2-40),[20,20]);
E0 = zeros(dim1,dim2);
M0 = zeros(dim1,dim2);

% creating gammas
gammas_ = Gammas(x1,N0,1,b_);

%     %plotting gammas
%         hold on
%         figure
%         surf(gammas_)
%         hold off

% creating vector of initial conditions
N0lin = reshape(N0,dim1*dim2,1);
E0lin = reshape(E0,dim1*dim2,1);
M0lin = reshape(M0,dim1*dim2,1);
y0 = [P0;N0lin;E0lin;M0lin];

% integrating diffeqs
[ts_vec,y_out] = ode45(@(t,y)ShSpSimdy(t,y,dim1,dim2),[0:0.1:days],y0);
n_ts = size(ts_vec,1);

% setting final values
P_out = y_out(:,1); % 2D matrices of cells per site (linear) per ts
N_outlin = y_out(:,2:dim1*dim2+1);
E_outlin = y_out(:,dim1*dim2+2:2*dim1*dim2+1);
M_outlin = y_out(:,2*dim1*dim2+2:end);

N_out = reshape(N_outlin,n_ts,dim1,dim2); % 3D matrices of cells per row per col per ts
E_out = reshape(E_outlin,n_ts,dim1,dim2);
M_out = reshape(M_outlin,n_ts,dim1,dim2);

N_tot = sum(N_outlin,2); % 1D vectors of total cells per ts
E_tot = sum(E_outlin,2);
M_tot = sum(M_outlin,2);

% plotting immune system response to pathogen
    figure
        semilogy(ts_vec,P_out,ts_vec,(M_tot + N_tot + E_tot))%,ts_vec,M_tot,ts_vec,N_tot,ts_vec,E_tot);
        axis([0 days 1 10^10])

% % plotting final N, M, and E populations
%         hold on
%         figure
%         surf(squeeze(N_out(n_ts,:,:)))
%         hold off
% 
%         hold on
%         figure
%         surf(squeeze(E_out(n_ts,:,:)))
%         hold off
% 
%         hold on
%         figure
%         surf(squeeze(M_out(n_ts,:,:)))
%         hold off
        
        
% % mutating infection loop
% 
% mRate = 
% 
% P02 = 10;
% N02 = N0 + M_out(n_ts,:,:);
% E02 = zeros(dim1,dim2);
% M02 = zeros(dim1,dim2);
% y02 = [P02;N02;E02;M02];
% 
% [ts_vec2,y_out2] = ode45(@(t,y)SimpleSimdy(t,y),[0:0.001:200],y02);
% n_ts2 = size(ts_vec2,1);
% 
% % setting final values
% P_out2 = y_out2(:,1);
% N_out2 = y_out2(:,2);
% E_out2 = y_out2(:,3);
% M_out2 = y_out2(:,4);
% 
% % plotting immune system response to pathogen
% figure
% semilogy(ts_vec2,P_out2,ts_vec2,(M_out2 + N_out2 + E_out2))%,ts_vec2,M_out2,ts_vec2,N_out2,ts_vec2,E_out2);
% axis([0 200 1 10^10])
%         



