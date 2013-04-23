% mutation 1D plotter

clear

global r_ h_ sigma_ de_ f_ c beta_ chi_ Qstep x0 dh_ ;
global b eps_ mu_ k_ Pdim1 Ldim1 Nstep Gamma_ ;

runnum = 2;
basecode = 'naive';
%datapath = 'C:\Documents and Settings\kimberly\Desktop\MATLAB\immune2012_data\'; %MOTHRA datapath
datapath = ['C:\Users\Kimberly\Google Drive\immunedata\' basecode '\'];%NEW laptop Gdrive
%datapath = 'C:\Users\Kimberly\dropbox\research\MATLAB\immune2012_data\';%laptop
%datapath = 'C:\Users\Kimberly\Desktop\Complex Systems\immune2012_data\'; %M-l transplant
afilename = [datapath 'a' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Nfilename = [datapath 'N' basecode num2str(runnum) '.txt'];
Efilename = [datapath 'E' basecode num2str(runnum) '.txt'];
Mfilename = [datapath 'M' basecode num2str(runnum) '.txt'];

% set parameters, read in days
params = setparams(afilename);
days = params{end,2};    % total days run & saved in file

% % gammas & lambdas (for beta_ landscape)
% b = 25;
% betas = [15;25;35];
% gammas1D = zeros(Pdim1,Ldim1);
% lambdas = zeros(Pdim1,size(betas,1));
% size(lambdas)
% for k=1:size(betas,1)
%     lambdas1D = zeros(Pdim1,1);
%     p_ = (1-exp(-1*((Pdim1)^2)/(8*(betas(k))^2)))^(-1);
%     for i=1:Pdim1;
%         lambdas1D(i) = 1 - p_*(1-exp(-1*((i-x0)^2)/(2*(betas(k))^2)));
%         for j=1:Ldim1;
%             gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
%         end
%     end
%     lambdas(:,k) = lambdas1D;
% end

% gammas & lambdas (for eps_ landscape)
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% % plot of affinity v. fitness curves (for beta_ landscape)
% figure
% hold on
% plot((1:Pdim1),gammas1D(:,x0+100),(1:Pdim1),lambdas(:,1),(1:Pdim1),lambdas(:,2),(1:Pdim1),lambdas(:,3))
% title('Affinity (\gamma(x,x_0)) v. Fitness (\lambda(x))')
% xlabel('position in shape space (site)')
% ylabel('value of affinity and fitness factors')
% legend(['\gamma(x,x_0 = 300), b = ' num2str(b)],['\lambda(x), \phi = ' num2str(betas(1))],['\lambda(x), \phi = ' num2str(betas(2))],['\lambda(x), \phi = ' num2str(betas(3))],'Location','Northwest')

% % plot of affinity v. fitness curve (for eps_ landscape)
% figure
% plot((1:Pdim1),gammas1D(:,x0),(1:Pdim1),lambdas1D)
% title(['Affinity (b = ' num2str(b) ') v. Fitness (\epsilon = ' num2str(eps_) ')'])
% xlabel('position in shape space (site)')
% ylabel('value of affinity and fitness factors')
% legend('Affinity \gamma(x,x_0)','Fitness \lambda(x)','Location','Northwest')
% 
% % plot of lambdas_ for several eps_ values
% epsvec_ = [1;5;15;25];
% lambdasvec_ = zeros(Pdim1,size(epsvec_,1));
% M = cell(size(epsvec_,1),1);
% for i=1:Pdim1;
%     for j=1:size(epsvec_,1)
%         lambdasvec_(i,j) = 1 - (2*epsvec_(j))/(Pdim1 + 2*epsvec_(j) - abs(Pdim1-2*i));
%     end
% end
% 
% figure
% for i=1:size(epsvec_,1)
%     size(lambdasvec_(i,:))
%     plot((1:Pdim1),lambdasvec_(:,i));
%     hold on;
%     hold all;
%     M{i} = ['\epsilon = ' num2str(Pdim1-epsvec_(i))];
% end
% title('Fitness Landscape \lambda(x)')
% xlabel('position in shape space (site)')
% %ylabel('value of affinity and fitness factors')
% ylabel('value of fitness function')
% legend(M{1},M{2},M{3},M{4},'Location','SouthEast')


% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Nplot = csvread(Nfilename);
Eplot = csvread(Efilename);
Mplot = csvread(Mfilename);

n_ts = size(tplot,1);

for i=1:size(Nplot,1)
    for j=1:size(Nplot,2)
        if Nplot(i,j) < mu_
            Nplot(i,j) = 0;
        end
        if Eplot(i,j) < mu_
            Eplot(i,j) = 0;
        end
        if Mplot(i,j) < mu_
            Mplot(i,j) = 0;
        end
    end
end


%plot of total pathogen v. total lymphocyte population
    Ptot = sum(Pplot,2);
    Ntot = sum(Nplot,2);
    Etot = sum(Eplot,2);
    Mtot = sum(Mplot,2);
    figure
    semilogy(tplot,Ptot,tplot,Ntot+Mtot+Etot,tplot,Ntot,tplot,Etot,tplot,Mtot)
    axis([0 days 1 10^10])
    title('Single-Infection Cell Populations, no mutation')%\phi = ' num2str(beta_)])
    title(['Single-Infection Cell Populations, \epsilon = ' num2str(eps_)])
    xlabel('duration of infection (days)')
    ylabel('total population (cells)')
    legend('Pathogen','All Lymphocytes','Naive','Effector','Memory','Location','NorthEast')
  
% plots of initial and final PNEM-distributions    
%     figure
%     plot((1:Pdim1),Pplot(229,:))
%    axis([0 Pdim1 0 12])
    
%     figure
%     plot((1:Pdim1),Pplot(end,:))
% %     
% %     figure 
%     plot((1:Pdim1),Nplot(end,:))
% 
%     figure 
%     plot((1:Pdim1),Eplot(end,:))
% 
%     figure 
%     plot((1:Pdim1),Mplot(end,:))
% 
% 
% plot of Psat over time (size = n_ts x Ldim1)
%     Pofy = zeros(n_ts,Ldim1);
%         for j = 1:Ldim1
%             Pofy(:,j)= sum(Pplot.*transpose(repmat(squeeze(gammas1D(:,j)),1,n_ts)),2);
%         end
%     Psat = Pofy./(k_.*ones(n_ts,Ldim1)+Pofy);
%     Xaxis = tplot;
%     Yaxis = (1:1:Ldim1);
%     figure
%     surf(Xaxis,Yaxis,transpose(Psat),'EdgeColor','none')
%     axis([0 days 0 Ldim1])
%     title('P_{sat} evolution over time')
%     ylabel('y-position in shape space')
%     xlabel('duration of infection (days)')
% 
% % 
% % contour plots of PNEM populations over time
% % NOTE these plots ARE ABSOLUTELY properly time-normalised
    Plog = ones(size(Pplot));
    Elog = ones(size(Eplot));
    for i=1:size(Plog,1)
        for j=1:size(Plog,2)
            if Pplot(i,j)>1
                Plog(i,j) = log(Pplot(i,j));
            end
            if Eplot(i,j)>1
                Elog(i,j) = log(Eplot(i,j));
            end
        end
    end
%     Xaxis = tplot;
%     Yaxis = (1:1:Pdim1);
%     figure
%     surf(Xaxis,Yaxis,transpose(Plog),'EdgeColor','none')
%     axis([0 days 0 Ldim1])
%     title(['Pathogen Evolution in Shape Space, b = ' num2str(b) ' (color on log scale)'])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
% %     v = [ mu_ 1 ];
% %     figure
% %     contourf(Xaxis,Yaxis,transpose(Pplot),v)
% %     axis([0 days 0 Pdim1])
% %     title(['Pathogen Evolution in Shape Space, \mu = ' num2str(mu_)])
% %     ylabel('position in shape space (site)')
% %     xlabel('duration of infection (days)')
% %     legend('Pathogen = \mu','Location','Northeast')
% % % 
% % %     Yaxis = (1:1:Ldim1);
% % %     figure
% % %     surf(Xaxis,Yaxis,transpose(Nplot),'EdgeColor','none')
% % %     axis([0 days 0 Ldim1])
% % %     title(['Naive Cell Evolution in Shape Space, \epsilon = ' num2str(eps_)])
% % %     ylabel('position in shape space (site)')
% % %     xlabel('duration of infection (days)')
% % %     size(Nplot)
% % %     Nplot(250,50)
% % %     
%     figure
%     surf(Xaxis,Yaxis,transpose(Nplot),'EdgeColor','none')
%     axis([0 days 0 Ldim1])
%     title(['Naive Cell Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     size(Nplot)
%     Nplot(250,50)
%     
%     figure
%     surf(Xaxis,Yaxis,transpose(Elog),'EdgeColor','none')
%     axis([0 days 0 400])
%     title(['Effector Evolution in Shape Space, \epsilon = ' num2str(eps_) ', color on log scale'])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     
%     figure
%     surf(Xaxis,Yaxis,transpose(Mplot),'EdgeColor','none')
%     axis([0 days 0 Ldim1])
%     title(['Memory Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     
%     figure
%     surf(Xaxis,Yaxis,transpose(Nplot+Mplot+Eplot),'EdgeColor','none')
%     axis([0 days 0 Ldim1])
%     title(['Total Lymphocyte Evolution in Shape Space, \epsilon = ' num2str(eps_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
% 
% plots of cutoff levels for whole infection
    v = [ mu_ 1 ];
%    v = [ 0 1 ];
    Xaxis = tplot;
    Yaxis = (1:1:Ldim1);
    figure
    contourf(Xaxis,Yaxis,transpose(Pplot),v)
    axis([0 days 0 Pdim1])
    title(['Pathogen Evolution in Shape Space, \mu = ' num2str(mu_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    legend('Pathogen = \mu','Location','Northeast')
    
    figure
    contourf(Xaxis,Yaxis,transpose(Eplot),v)
    axis([0 days 0 Ldim1])
    title(['Effector Evolution in Shape Space, \mu = ' num2str(mu_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    legend('Effector = \mu','Location','Northeast')
    
    figure
    contourf(Xaxis,Yaxis,transpose(Mplot),v)
    axis([0 days 0 Ldim1])
    title(['Memory Evolution in Shape Space, \mu = ' num2str(mu_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    legend('Memory = \mu','Location','Northeast')

    figure
    contourf(Xaxis,Yaxis,transpose(Nplot),v)
    axis([0 days 0 Ldim1])
    title(['Naive Evolution in Shape Space, \mu = ' num2str(mu_)])
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    legend('Naive = \mu','Location','Northeast')

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
