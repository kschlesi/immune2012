% mutation 1D plotter

clear

% global r_ h_ sigma_ de_ f_ k_ c ;
global b beta_ mu_;

%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\';
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';
tfilename = [datapath 'tbound17.txt'];
Pfilename = [datapath 'Pbound17.txt'];
Nfilename = [datapath 'Nbound17.txt'];
Efilename = [datapath 'Ebound17.txt'];
Mfilename = [datapath 'Mbound17.txt'];

days = 1000;       % total days run

% dimensions of 1D shape space
Pdim1 = 400;
Ldim1 = 400;
x0 = 200;

% % gammas & lambdas
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
p_ = (1-exp(-1*((Pdim1)^2)/(8*beta_^2)))^(-1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - p_*(1-exp(-1*((i-x0)^2)/(2*beta_^2)));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% figure
% plot((1:Pdim1),gammas1D(:,x0),(1:Pdim1),lambdas1D)
% title(['Affinity (b = ' num2str(b) ') v. Fitness (\phi = ' num2str(beta_) ')'])
% xlabel('position in shape space (site)')
% ylabel('value of affinity and fitness factors')
% legend('Affinity \gamma(x,x_0)','Fitness \lambda(x)','Location','Northwest')





% % data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(tplot);


% % plot of total pathogen v. total lymphocyte population
    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    figure
    semilogy(tplot,Ptot,tplot,Ntot+Mtot+Etot)
    axis([0 days 1 10^9])
    title(['Single-Infection Cell Populations, \phi = ' num2str(beta_)])
    xlabel('duration of infection (days)')
    ylabel('total population (cells)')
    legend('Pathogen','All Lymphocytes','Location','Northeast')

    
% % % plot of initial and final P-distributions    
% %     figure
% %     plot((1:Pdim1),Pplot(1,:))
% %     axis([0 Pdim1 0 12])
% %     
% %     figure
% %     plot((1:Pdim1),Pplot(end,:))
% %     
% %     figure 
% %     plot((1:Pdim1),Nplot(end,:))
% % 
% %     figure 
% %     plot((1:Pdim1),Eplot(end,:))
% % 
% %     figure 
% %     plot((1:Pdim1),Mplot(end,:))
% %     
% % contour plots of PNEM populations over time
% % NOTE these plots ARE ABSOLUTELY properly time-normalised
    
    Yaxis = tplot;
    Xaxis = (1:1:Pdim1);
    figure
    v = [1 10:50000:100000000];
    contour(Xaxis,Yaxis,Pplot,v)
    axis([0 400 0 days])
    title(['Pathogen Evolution in Shape Space, \phi = ' num2str(beta_)])
    xlabel('position in shape space (site)')
    ylabel('duration of infection (days)')
%     v = [ mu_ 1 ];
%     figure
%     contour(Xaxis,Yaxis,Pplot,v)
%     title(['Pathogen Evolution in Shape Space, \phi = ' num2str(beta_)])
%     xlabel('position in shape space (site)')
%     ylabel('duration of infection (days)')
%     legend('Pathogen = \mu','Location','Northeast')

% %     Xaxis = (1:1:Ldim1);
% %     figure
% %     v = [0:0.5:3];
% %     contour(Xaxis,Yaxis,Nplot,v)
% %     
% %     figure
% %     v = [1 10 50 100 200 300 500:500:10000];
% %     contour(Xaxis,Yaxis,Eplot,v)
% % 
% %     figure
% %     v = [1 10 50 100 200 300 500:500:10000];
% %     contour(Xaxis,Yaxis,Mplot,v)
