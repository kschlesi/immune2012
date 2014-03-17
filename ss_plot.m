% mutation 1D plotter

clear all -except mesh_tests

runnum = 917;
basecode = 'cmeshtwo';
datapath = ['/Users/kimberly/Google Drive/immunedata/PL13/' basecode '/'];
bfilename = [datapath 'b' basecode num2str(runnum) '.txt'];
tfilename = [datapath 't' basecode num2str(runnum) '.txt'];
Pfilename = [datapath 'P' basecode num2str(runnum) '.txt'];
Lfilename = [datapath 'L' basecode num2str(runnum) '.txt'];

% set parameters
days = 0;
params = setparams(bfilename);
for i=1:size(params,1)
    eval([char(params{i,1}) ' = ' num2str(params{i,2}) ';']);
    eval([char(params{i,1}) 'units = char(params{i,3});']);
end
clear params;

% data and time vector
tplot = csvread(tfilename);
Pplot = csvread(Pfilename);
Lplot = csvread(Lfilename);

n_ts = size(tplot,1);

% update paramfile with last 'days' value saved if ode45 was interrupted
if ~days
    tend = cell(1,3);                   
    tend{1,1} = 'days';                
    tend{1,2} = tplot(end);
    tend{1,3} = 'days';
    cell2csv(bfilename,tend,1); 
    clear tend;
    days = tplot(end);
end    
    
Pplot = Pplot.*(Pplot>=mu_);
Lplot = Lplot.*(Lplot>=mu_);
Ptot = sum(Pplot,2);
Ltot = sum(Lplot,2);

% %%%%%%%%%%%% plots of ratios to determine chronicity....
% gammas1D = zeros(Pdim1,Ldim1);   % matrix of affinities
% lambdas1D = Lambdas(eps_,Pdim1); % vector of pathogen fitnesses      
% for i=1:Pdim1;
%     gammas1D(i,:) = Gammas([i,1],ones(Ldim1,1),1,b);
% end
% 
% avgPx = Pplot*(1:400)'./Ptot;
% figure
% plot(tplot,avgPx)
% gammas_at_avgPx = zeros(size(Lplot));
% for i=1:size(tplot,1)
%     gammas_at_avgPx(i,:) = gammas1D(round(avgPx(i)),:);
% end
% omegas_at_avgPx = sum(gammas_at_avgPx.*Lplot,2);
% omegasPtot_number = zeros(size(omegas_at_avgPx));
% for time=1:size(omegas_at_avgPx,1)
%     omegasPtot_number = lambdas1D(round(avgPx(time)))*(1-Ptot(time)/K_)/h_;
% end
% figure
% plot(tplot,omegas_at_avgPx./omegasPtot_number)

% calculate maxd from first Ppeak; compare with estimate
Pderiv = diff(Ptot);
peaks = (Pderiv<0).*(circshift(Pderiv,1)>=0);
Ppeak = Ptot(find(peaks,1,'first'));
%disp(Ptot(1:end-1).*peaks);
%peakind = (derivtest<0).*(circshift(derivtest,1)>0);
if Ppeak
peakindex = find(Ptot==Ppeak);
peaktime = tplot(peakindex);
maxsite = find(Pplot(peakindex+2,:)==0,1,'first')-1;
maxd = maxsite-x0;
curlyL = sqrt(2/pi)/(Pdim1*chi_);
disp([maxd sqrt(curlyL*Ppeak*(1-Ppeak/K_))]);
end

% Pstrains = sum(Pplot>0,2);
% figure
% plot(tplot,Pstrains,'b--',tplot,[0;diff(Pstrains)],'r')
% axis([0 days 0 Pdim1])
%axis([50 500 -5 25])


% day1 = 600;
% day2 = 800;
% indx = find(tplot>day1);
% indx1 = indx(1)-1;
% indx = find(tplot>day2);
% indx2 = indx(1)-1;
% dispel=[sum(Lplot(indx1,:)) sum(Lplot(indx2,:))];
% disp(dispel);
% disp(dispel(2)/dispel(1));
    
%     figure
%     plot((1:1:size(tplot,1)),tplot)
%     
%     Ptot(end)
%     Ltot(end)
%   
%     figure
%     plot((1:1:Ldim1),Lplot(1500,:))
%     axis([50 100 0 120000])
    
% plots of initial and final PL-distributions    
    
%     figure
%     plot((1:Pdim1),Pplot(end,:))

%     figure 
%     plot((1:Ldim1),Lplot(1,:))
% 
%     figure 
%     plot((1:Ldim1),Lplot(end,:))
% 
%     Ltot(1)/Ldim1
%     Ltot(end)/Ldim1
%     Gamma_/delta_
 
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
% % % contour plots of PL populations over time
% % % NOTE these plots ARE ABSOLUTELY properly time-normalised
    Xaxis = tplot;
    Yaxis = (1:1:Pdim1);
    logsurf(Xaxis,Yaxis,Pplot')
%     contdata = Pplot;
%     contdata(Pplot<mu_) = -Inf;
%     figure
%     v = [mu_,500,1e07:1e07:10e07];
%     vv = [ mu_ 1 ];
%     contourf(Xaxis,Yaxis,contdata',vv)
%     hold on
%     contourf(Xaxis,Yaxis, contdata',v,'EdgeColor','none')
%     clabel(C)
    set(gca,'GridLineStyle','none')
    axis([0 days 0 Pdim1])
    axis([50 2000 0 300])
    %title('Pathogen Evolution in Shape Space')
    %ylabel('position in shape space (site)')
    %xlabel('duration of infection (days)')
%     %clear Pplot;
    
%     %plot of total pathogen v. total lymphocyte population
    figure
    semilogy(tplot,Ptot,'-',tplot,Ltot,'-')
    %hold on
    %plot(tplot,(Ltot(1)).*ones(size(tplot)),'r')
    %plot(tplot,(R_).*ones(size(tplot)),'--r')
    axis([0 days 1 10^10])
    axis([0 25 1 10^10])
    title('Single-Infection Cell Populations')
%    title(['Single-Infection Cell Populations, b = ' num2str(b)])
    xlabel('duration of infection (days)')
    ylabel('total population (cells)')
    legend('Pathogen','Lymphocytes','Location','NorthWest')
    %clear Ptot;
    %clear Ltot;
%     
%     figure
%     semilogy(tplot,Pplot(:,6))
%     hold on
%     hold all
%     semilogy(tplot,Lplot(:,6))
%     title([basecode num2str(runnum)])
%     %axis([0 10 0 2.6e6])

    
%     Plog = ones(size(Pplot));
%     Llog = ones(size(Lplot));
%     for i=1:size(Plog,1)  % log-scaling the values
%         for j=1:size(Plog,2)
%             if Pplot(i,j)>1
%                 Plog(i,j) = log(Pplot(i,j));
%             else
%                 Plog(i,j) = 0;
%             end
%             if Lplot(i,j)>1
%                 Llog(i,j) = log(Lplot(i,j));
%             else
%                 Plog(i,j) = 0;
%             end
%         end
%     end
%     
%     Xaxis = tplot;
%     Yaxis = (1:1:Pdim1);
%     figure
%     surf(Xaxis,Yaxis,transpose(Plog),'EdgeColor','none')
%     axis([0 days 0 Pdim1])
%     %axis([0 20 0 100])
%     title(['Pathogen Evolution in Shape Space, with mutation' '(color on log scale)'])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     colorbar('Location','EastOutside')

    Yaxis = (1:1:Ldim1);
    logsurf(Xaxis,Yaxis,Lplot')
    axis([0 days 0 Ldim1])
    %axis([0 10 0 100])
    title('Lymphocyte Evolution in Shape Space')
    ylabel('position in shape space (site)')
    xlabel('duration of infection (days)')
    %clear Lplot;
      
% % plots of cutoff levels for whole infection
%     v = [ mu_ 1 ];
%     Xaxis = tplot;
%     Yaxis = (1:1:Ldim1);
%     figure
%     contourf(Xaxis,Yaxis,Pplot',v)
%     axis([0 days 0 Pdim1])
%     title(['Pathogen Evolution in Shape Space, \mu = ' num2str(mu_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     legend('Pathogen = \mu','Location','Northeast')
%     
%     figure
%     contourf(Xaxis,Yaxis,transpose(Lplot),v)
%     axis([0 days 0 Ldim1])
%     title(['Lymphocyte Evolution in Shape Space, \mu = ' num2str(mu_)])
%     ylabel('position in shape space (site)')
%     xlabel('duration of infection (days)')
%     legend('Lymphocytes = \mu','Location','Northeast')

ix=100;
Pline = Pplot(find(tplot>ix,1,'first'),:);
Lline = Lplot(find(tplot>ix,1,'first'),:);

figure    % plot of P0 and L0 distributions at days+olddays
hold on
hold all
%plot((1:1:Pdim1),Pplot(end,:))
%plot((1:1:Ldim1),Lplot(end,:))
%plot((1:1:Ldim1),Lplot(1,:))
%plot((1:1:Ldim1),R_/Ldim1)
%title(['P0 and L0 distributions at t = ' num2str(days) ' days'])
if Pline
    if Lline
plot((1:1:Pdim1),Pline(end,:))
plot((1:1:Ldim1),Lline(end,:))
    end
end
title(['P0 and L0 distributions at t = ' num2str(ix) ' days'])
xlabel('location in shape space (site)')
ylabel('population (cells/\mul)')
legend('Pathogen','Lymphocytes')
%axis([0 100 0 2.4e4])
    
%%%%%%%%%%%%%%gammas, lambdas, etc
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

% % gammas & lambdas (for eps_ landscape)
% gammas1D = zeros(Pdim1,Ldim1);
% lambdas1D = zeros(Pdim1,1);
% for i=1:Pdim1;
%     lambdas1D(i) = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
%     for j=1:Ldim1;
%         gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
%     end
% end

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
% epsvec_ = [1;4;15;25];
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
%     M{i} = ['\epsilon = ' num2str((Pdim1-2*epsvec_(i))/Pdim1)];
% end
% axis([1 Pdim1 0 1.1])
% title('Fitness Landscape \lambda(x)')
% xlabel('position in shape space (site)')
% %ylabel('value of affinity and fitness factors')
% ylabel('value of fitness function')
% legend(M{1},M{2},M{3},M{4},'Location','SouthEast')

