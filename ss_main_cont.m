% mutation1dall_cont

clear

global r_ sigma_ de_ f_ k_ c gammas1D h_ ;

r_ = 3.3;
h_ = 10^-3;
sigma_ = 3;
de_ = 0.35;
k_ = 10^5;
f_ = 0.1;
c = 0.05;
days = 5;
stepsize = 0.1;
olddays = 10;
oldss = 0.1;
%mrate = 0.7; % per cell per day
b = 15;
N0density = 3;
x0 = 500;

oldt_vec = (0:oldss:olddays);
old_ts = size(oldt_vec,2);

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


% initial inoculation in shape space
P0in = csvread('Pdif10a.txt');
N0in = csvread('Ndif10a.txt');
E0in = csvread('Edif10a.txt');
M0in = csvread('Mdif10a.txt');
P0 = shiftdim(P0in(old_ts,:),1);
N0 = shiftdim(N0in(old_ts,:),1);
E0 = shiftdim(E0in(old_ts,:),1);
M0 = shiftdim(M0in(old_ts,:),1);


% creating initial conditions vector

y0 = [P0;N0;E0;M0];

figure
 semilogy((1:Pdim1+3*Ldim1),(shiftdim(y0)))
% semilogy((1:Ldim1),(shiftdim(N0)))
% semilogy((1:Ldim1),(shiftdim(E0)))
% semilogy((1:Ldim1),(shiftdim(M0)))


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
dlmwrite('Pnewbie.txt',P_out);
concat('Pdif10a.txt','Pnewbie.txt');
dlmwrite('Nnewbie.txt',N_out);
concat('Ndif10a.txt','Nnewbie.txt');
dlmwrite('Enewbie.txt',E_out);
concat('Edif10a.txt','Enewbie.txt');
dlmwrite('Mnewbie.txt',M_out);
concat('Mdif10a.txt','Mnewbie.txt');


% plot initial & final distributions
    figure
    plot((1:1:Pdim1),P0)
    
    figure
    Pfin = squeeze(P_out(n_ts,:));
    plot((1:1:Pdim1),Pfin)
    


    Ptot = sum(P_out,2);
    Ntot = sum(N_out,2);
    Etot = sum(E_out,2);
    Mtot = sum(M_out,2);
    figure
    semilogy(ts_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
    axis([0 days 1 10^8])
    
    
        figure
    hold on
    surf(P_out,'MeshStyle','row')
    hold off
    axis([0 Pdim1 0 n_ts 0 max(P0(:,x0))])
%    set(gca,'ZScale','log')

    figure
    hold on
    surf(N_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 max(N_out(:,x0))])

    figure
    hold on
    surf(E_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 max(E_out(:,x0))])
    
    figure
    hold on
    surf(M_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 max(M_out(:,x0))])
