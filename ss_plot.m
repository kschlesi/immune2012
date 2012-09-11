% mutation 1D plotter

clear

global r_ h_ sigma_ de_ f_ k_ c b p_ beta_ ;

Pfilename = 'Ploop3.txt';
Nfilename = 'Nloop3.txt';
Efilename = 'Eloop3.txt';
Mfilename = 'Mloop3.txt';

days = 1;        % total days run
stepsize = 0.1;   % interval (days) at which ode45 was called

% dimensions of 1D shape space
Pdim1 = 600;
Ldim1 = 600;
x0 = 300;

% gammas & lambdas
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = p_*(1-exp(-1*((i-x0)^2)/(2*beta_^2)));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% data and time vector
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(Pplot,1);
ts_vec = (0:n_ts-1);            % vector of all internal timesteps taken by ode45
plot_vec = ts_vec.*days./n_ts;  % ts_vec rescaled to units of days (for plotting)


% plot of total pathogen v. total lymphocyte population
    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    figure
    semilogy(plot_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
    axis([0 days 1 10^8])
    

% contour plot of P population over time
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