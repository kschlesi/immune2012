% movie of given run(s) in sssDiffuseDy

clear

global b_ r_ h_ sigma_ k_ de_ f_ gamma_lib;

% set attributes of run to be filmed
days = 20;
stepsize = 0.1;
ts_v = [0:stepsize:days];
timesteps = size(ts_v,2); %number of timesteps in run

Pdim1 = 60;
Pdim2 = 48;
Ldim1 = 60;
Ldim2 = 48;

% read in & array data (timestep x dim1 x dim2)
Plin = csvread('Pout14.txt');
Nlin = csvread('Nout14.txt');
Elin = csvread('Eout14.txt');
Mlin = csvread('Mout14.txt');
Pplot = reshape(Plin,timesteps,Pdim1,Pdim2);
Nplot = reshape(Nlin,timesteps,Ldim1,Ldim2);
Eplot = reshape(Elin,timesteps,Ldim1,Ldim2);
Mplot = reshape(Mlin,timesteps,Ldim1,Ldim2);

% figure
% surf(squeeze(Eplot(timesteps,:,:)))
% set(gca,'ZScale','log')
% axis([0 48 0 60 1 4e4 1 400])
% get(gca,'ZScale')


% sean's matlab movie code
figure
surf(squeeze(Pplot(1,:,:)));
set(gca,'ZScale','log')
axis([0 48 0 60 1 10^10 1 100])
set(gca,'NextPlot','replacechildren','ZScale','log');
% Preallocate the struct array for the struct returned by getframe
F(20) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:timesteps
    surf(squeeze(Pplot(j,:,:)))
    hold on
    set(gca,'ZScale','log')
    axis([0, 48, 0, 60, 1, 10^10, 1, 100])
    hold off
    F(j) = getframe;
end

figure
surf(squeeze(Nplot(1,:,:)));
axis([0 48 0 60 1 4 1 100])
set(gca,'NextPlot','replacechildren');
% Preallocate the struct array for the struct returned by getframe
G(20) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:timesteps
    surf(squeeze(Nplot(j,:,:)))
    hold on
    axis([0, 48, 0, 60, 1, 4, 1, 100])
    hold off
    G(j) = getframe;
end

figure
surf(squeeze(Eplot(1,:,:)));
set(gca,'ZScale','log')
axis([0 48 0 60 1 10^5 1 100])
set(gca,'NextPlot','replacechildren','ZScale','log');
% Preallocate the struct array for the struct returned by getframe
H(20) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:timesteps
    surf(squeeze(Eplot(j,:,:)))
    hold on
    set(gca,'ZScale','log')
    axis([0, 48, 0, 60, 1, 10^5, 1, 100])
    hold off
    H(j) = getframe;
end

figure
surf(squeeze(Mplot(1,:,:)));
axis([0 48 0 60 1 4000 1 100])
set(gca,'NextPlot','replacechildren');
% Preallocate the struct array for the struct returned by getframe
I(20) = struct('cdata',[],'colormap',[]);
% Record the movie
for j = 1:timesteps
    surf(squeeze(Mplot(j,:,:)))
    hold on
    axis([0, 48, 0, 60, 1, 4000, 1, 100])
    hold off
    I(j) = getframe;
end


movie(F)
movie(G)
movie(H)
movie(I)