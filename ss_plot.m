% mutation 1D plotter

clear

% global r_ h_ sigma_ de_ f_ c ;
global b eps_ mu_ k_ ;

runnum = 9;
basecode = 'quant';
datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';%laptop
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %M-l transplant
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

days = 10;       % total days run

% dimensions of 1D shape space
Pdim1 = 400;
Ldim1 = 400;
x0 = 200;

% gammas & lambdas
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% figure
% plot((1:Pdim1),gammas1D(:,x0),(1:Pdim1),lambdas1D)
% title(['Affinity (b = ' num2str(b) ') v. Fitness (\epsilon = ' num2str(eps_) ')'])
% xlabel('position in shape space (site)')
% ylabel('value of affinity and fitness factors')
% legend('Affinity \gamma(x,x_0)','Fitness \lambda(x)','Location','Northwest')





% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(tplot,1);


% plot of total pathogen v. total lymphocyte population
    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    figure
    semilogy(tplot,Ptot,tplot,Ntot+Mtot+Etot)
    axis([0 days 1 10^10])
    title(['Single-Infection Cell Populations, \epsilon = ' num2str(eps_)])
    xlabel('duration of infection (days)')
    ylabel('total population (cells)')
    legend('Pathogen','All Lymphocytes','Location','NorthWest')

    
% % plot of initial and final PNEM-distributions    
%     figure
%     plot((1:Pdim1),Pplot(1,:))
%     axis([0 Pdim1 0 12])
%     
%     figure
%     plot((1:Pdim1),Pplot(end,:))
%     
%     figure 
%     plot((1:Pdim1),Nplot(end,:))
% 
%     figure 
%     plot((1:Pdim1),Eplot(end,:))
% 
%     figure 
%     plot((1:Pdim1),Mplot(end,:))


% plot of Psat over time (size = n_ts x Ldim1)
    Pofy = zeros(n_ts,Ldim1);
        for j = 1:Ldim1
            Pofy(:,j)= sum(Pplot.*transpose(repmat(squeeze(gammas1D(:,j)),1,n_ts)),2);
        end
    Psat = Pofy./(k_.*ones(n_ts,Ldim1)+Pofy);
    Xaxis = tplot;
    Yaxis = (1:1:Ldim1);
    v = (0:0.05:1);
    figure
    contourf(Xaxis,Yaxis,transpose(Psat),v)
    axis([0 days 0 Ldim1])
    title('P_{sat} evolution over time')
    ylabel('y-position in shape space')
    xlabel('duration of infection (days)')
    legend


% contour plots of PNEM populations over time
% NOTE these plots ARE ABSOLUTELY properly time-normalised
    Xaxis = tplot;
    Yaxis = (1:1:Pdim1);
    figure
    v = [1 10:50000:10000000];
    contourf(Xaxis,Yaxis,transpose(Pplot),v)
    axis([0 days 0 Pdim1])
    title(['Pathogen Evolution in Shape Space, \epsilon = ' num2str(eps_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    v = [ mu_ 1 ];
    figure
    contourf(Xaxis,Yaxis,transpose(Pplot),v)
    axis([0 days 0 Pdim1])
    title(['Pathogen Evolution in Shape Space, \epsilon = ' num2str(eps_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    legend('Pathogen = \mu','Location','Northeast')

    Yaxis = (1:1:Ldim1);
    figure
%    v = (0:0.5:3);
    contourf(Xaxis,Yaxis,transpose(Nplot),25)
    axis([0 days 0 Ldim1])
    title(['Naive Cell Evolution in Shape Space, \epsilon = ' num2str(eps_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    
    figure
%    v = [0 1 10 50 100 200 300 500:500:10000];
    contourf(Xaxis,Yaxis,transpose(Eplot),25)
    axis([0 days 0 Ldim1])
    title(['Effector Evolution in Shape Space, \epsilon = ' num2str(eps_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    
    figure
    v = [0 1 10 50 100 200 300 500:500:10000];
    contourf(Xaxis,Yaxis,transpose(Mplot),25)
    axis([0 days 0 Ldim1])
    title(['Memory Evolution in Shape Space, \epsilon = ' num2str(eps_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')

    
% contourf plots, of PNEM evolution over time, normalised by total number
% of cells at each timestep

%     Ptotal = repmat(Ptot,1,Pdim1);
%     Ntotal = repmat(Ntot,1,Ldim1);
%     Etotal = repmat(Etot,1,Ldim1);
%     Mtotal = repmat(Mtot,1,Ldim1);
% 
%     Xaxis = tplot;
%     Yaxis = (1:1:Pdim1);
%     figure
% %    v = [0:0.01:1];
%     contourf(Xaxis,Yaxis,transpose(Pplot./Ptotal),20)
%     axis([0 days 0 Pdim1])
%     title(['Normalised Pathogen Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
% 
%     Yaxis = (1:1:Ldim1);
%     figure
% %    v = (0:0.5:3)/3;
%     contourf(Xaxis,Yaxis,transpose(Nplot./Ntotal),20)
%     axis([0 days 0 Ldim1])
%     title(['Normalised Naive Cell Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     
%     figure
%     contourf(Xaxis,Yaxis,transpose(Eplot./Etotal),20)
%     axis([0 days 0 Ldim1])
%     title(['Normalised Effector Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     
%     figure
%     contourf(Xaxis,Yaxis,transpose(Mplot./Mtotal),20)
%     axis([0 days 0 Ldim1])
%     title(['Normalised Memory Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
