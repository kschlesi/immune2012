% mutation 1D plotter

clear

r_ = 3.3;
h_ = 10^-5;
sigma_ = 3;
de_ = 0.35;
k_ = 10^5;
f_ = 0.1;
D = 10;
days = 12;
stepsize = 0.1;
%mrate = 0.7; % per cell per day
b = 50;
N0density = 3;
x0 = 100;

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

% data and time vector
Pplot = csvread('Pdif9c.txt');
Nplot = csvread('Ndif9c.txt');
Eplot = csvread('Edif9c.txt');
Mplot = csvread('Mdif9c.txt');

ts_vec = (0:stepsize:days);
n_ts = size(ts_vec,2);

% figure
% hold on
% surf(Pplot)
% hold off




    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    figure
    semilogy(ts_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
    axis([0 days 1 10^8])
    


    figure
    col(:,:,1) = rand(n_ts,Pdim1);
    col(:,:,2) = rand(n_ts,Pdim1);
    col(:,:,3) = rand(n_ts,Pdim1);
    v = [1 10 50 100 200 300 500:500:10000];
    contour(Pplot,v)%,'MeshStyle','row')
%    axis([0 Pdim1 0 n_ts 0 10^5])
%    %set(gca,'CLim',[30 85690]);
%    v=caxis;
%    v
%    set(gca,'ZScale','log')

%     figure
%     hold on
%     surf(Nplot)%,'MeshStyle','row')
%     hold off
%     axis([0 Ldim1 0 n_ts 0 N0density])
% 
%     figure
%     hold on
%     surf(Eplot)%,'MeshStyle','row')
%     hold off
%     axis([0 Ldim1 0 n_ts 0 10^4])
%     
%     figure
%     hold on
%     surf(Mplot)%,'MeshStyle','row')
%     hold off
%     axis([0 Ldim1 0 n_ts 0 10^3])
% 
% 
