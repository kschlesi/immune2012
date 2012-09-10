% for reading and plotting results of simulations 
% saved in arrays from sssim_diffuse.m

clear

global b_ r_ h_ sigma_ k_ de_ f_ gamma_lib;

days = 20;
stepsize = 0.1;
ts_v = [0:stepsize:days];
timesteps = size(ts_v,2); %number of timesteps in run to be read

Pdim1 = 60;
Pdim2 = 48;
Ldim1 = 60;
Ldim2 = 48;

% read in and reshape saved vectors
Plin = csvread('Pout17.txt');
Nlin = csvread('Nout17.txt');
Elin = csvread('Eout17.txt');
Mlin = csvread('Mout17.txt');

Pplot = reshape(Plin,timesteps,Pdim1,Pdim2);
Nplot = reshape(Nlin,timesteps,Ldim1,Ldim2);
Eplot = reshape(Elin,timesteps,Ldim1,Ldim2);
Mplot = reshape(Mlin,timesteps,Ldim1,Ldim2);

% plot, referencing data by Xplot(timestep,ssloc_x,ssloc_y)

    %satfunc(Peff) at particular y-points
     repgammaliby = shiftdim(repmat(squeeze(gamma_lib(:,:,30,25)),[1,1,timesteps]),2);
     Pofy = squeeze(sum(sum(Pplot.*repgammaliby,2),3))./(Pdim1*Pdim2);
     newsatfunc1 = Pofy./(k_*ones(timesteps,1)+Pofy);
%     repgammaliby = shiftdim(repmat(squeeze(gamma_lib(:,:,45,20)),[1,1,timesteps]),2);
%     Pofy = squeeze(sum(sum(Pplot.*repgammaliby,2),3))./(Pdim1*Pdim2);
%     newsatfunc2 = Pofy./(k_*ones(timesteps,1)+Pofy);
%      figure
%      semilogy(ts_v,newsatfunc1)
%     

        
    % total P, N, E, M over time
    Ptot = squeeze(sum(sum(Pplot,2),3));
    Ntot = squeeze(sum(sum(Nplot,2),3));
    Etot = squeeze(sum(sum(Eplot,2),3));
    Mtot = squeeze(sum(sum(Mplot,2),3));
    Ltot = Ntot+Etot+Mtot;
      figure
      semilogy(ts_v,Ptot,ts_v,Ltot)%,ts_v,Ntot,ts_v,Etot,ts_v,Mtot)
      axis([0 days 1 10^10])
    
   %initial P     
         hold on
         figure
         surf(squeeze(Pplot(1,:,:)))
         hold off

   %final N, E, M
         hold on
         figure
         surf(squeeze(Nplot(timesteps,:,:)))
         hold off
         
         hold on
         figure
         surf(squeeze(Eplot(timesteps,:,:)))
         hold off
         
         hold on
         figure
         surf(squeeze(Pplot(timesteps,:,:)))
         hold off
