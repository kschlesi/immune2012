% movie of given run(s) in P-only mutation script

clear

global b_ r_ h_ sigma_ k_ de_ f_ gamma_lib;

% set attributes of run to be filmed
days = 5;
stepsize = 0.1;
ts_v = [0:stepsize:days];
timesteps = size(ts_v,2); %number of timesteps in run

Pdim1 = 20;
Pdim2 = 20;

% read in & array data (timestep x dim1 x dim2)
Plin = csvread('Pmut5c.txt');
Pplot = reshape(Plin,timesteps,Pdim1,Pdim2);

Ptot = squeeze(sum(sum(Pplot,3),2));
figure
semilogy(ts_v,Ptot)

% figure
% surf(squeeze(Eplot(timesteps,:,:)))
% set(gca,'ZScale','log')
% axis([0 48 0 60 1 4e4 1 400])
% get(gca,'ZScale')


% sean's matlab movie code
figure
surf(squeeze(Pplot(1,:,:)));
%set(gca,'ZScale','log')
axis([0 Pdim2 0 Pdim1 0 10^8])
caxis auto
set(gca,'NextPlot','replacechildren');
% Preallocate the struct array for the struct returned by getframe
F(3) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:timesteps
    surf(squeeze(Pplot(j,:,:)))
    hold on
%    set(gca,'ZScale','log')
    axis([0, Pdim2, 0, Pdim1, 0, 10^8])
    caxis auto
    hold off
    F(j) = getframe;
end
movie(F,5)
