% test for diffused antigen mutation in shape space

clear

global r_ days stepsize mrate ijs;

r_ = 1;
days = 5;
stepsize = 0.1;
mrate = 0.7; % per cell per day

% dimensions of 2D shape space
Pdim1 = 20; 
Pdim2 = 20;

% center and max amount of initial gaussian inoculation in shape space
x0 = [10,10];
Pmax0 = 10;
Pdiff0 = 2;

% setting initial conditions for P
P0 = Gammas(x0,zeros(Pdim1,Pdim2),Pmax0,Pdiff0);  % initial gaussian distribution of pathogen
y0 = reshape(P0,Pdim1*Pdim2,1);

ijs = zeros(Pdim1,Pdim2,2); % every Pdim site holds ordered pair of its own location
for i=1:Pdim1
    for j=1:Pdim2
        ijs(i,j,:) = [i j];
    end
end
ijs = cast(ijs,'int16');

% integrating P diffeq in time
options = odeset('AbsTol',1e-3);
[ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Pdim2),[0:stepsize:days],y0,options);
n_ts = size(ts_vec,1);
P_out = reshape(y_out,n_ts,Pdim1,Pdim2);

% save results!!!
dlmwrite('Pmut4c.txt',y_out);

% plot initial & final P distributions
    hold on
    figure
    surf(P0)
    hold off

    figure
    hold on
    Pfin = squeeze(P_out(n_ts,:,:));
    surf(Pfin)
    hold off
