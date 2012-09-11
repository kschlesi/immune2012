% 1D fitness/mutation test
% with cutoff implemented by looping calls to ode45

clear

global r_ h_ sigma_ de_ f_ k_ c p_ beta_ lambdas1D gammas1D ;

days = 1;
stepsize = 0.1;

Pfilename = 'Ploop2.txt';
Nfilename = 'Nloop2.txt';
Efilename = 'Eloop2.txt';
Mfilename = 'Mloop2.txt';

% setting necessary parameters
r_ = 3.3;
h_ = 10^-3;
sigma_ = 3;
de_ = 0.35;
k_ = 10^5;
f_ = 0.1;
c = 0.5;
b = 25;
beta_ = 40; 


% dimensions of 1D shape space
Pdim1 = 600;
Ldim1 = 600;
x0 = 300;

% gammas & lambdas
p_ = (1-exp(-1*((Pdim1)^2)/(8*beta_^2)));
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = p_*(1-exp(-1*((i-x0)^2)/(2*beta_^2)));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% initial configurations
Pmax0 = 10;
Pdiff0 = 12;
G0 = Pmax0.*ones(Pdiff0/2,1);
P0 = padarray(G0,Pdim1/2-Pdiff0/4,'both');

N0density = 3;
N0 = N0density.*ones(Ldim1,1);
E0 = zeros(Ldim1,1);
M0 = zeros(Ldim1,1);

% P0 = zeros(Pdim1,1); % initial gaussian distribution of pathogen
% for i=1:Pdim1;
%     P0(i) = Pmax0*exp(-1*((i-x0)^2)/(2*Pdiff0^2));
% end

% creating initial conditions vector
y0 = [P0;N0;E0;M0];

dlmwrite(Pfilename,transpose(P0));
dlmwrite(Nfilename,transpose(N0));
dlmwrite(Efilename,transpose(E0));
dlmwrite(Mfilename,transpose(M0));    


% integrating diffeqs in time with a FOR LOOP
options = odeset('AbsTol',1e-3);
nsteps = cast(days/stepsize,'uint16');
n_ts = 0;
for j=1:nsteps
    
    % integrate between two external steps (of size stepsize)
    i = cast(j,'double');
    [ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),[(i-1)*stepsize,i*stepsize],y0,options);

    % add new internal steps to overall n_ts
    n_ts = n_ts + size(ts_vec,1)-1;
    
    % save & append y-output
    P_out = y_out(:,1:Pdim1);
    N_out = y_out(:,Pdim1+1:Pdim1+Ldim1);
    E_out = y_out(:,Pdim1+Ldim1+1:Pdim1+2*Ldim1);
    M_out = y_out(:,Pdim1+2*Ldim1+1:end);

    dlmwrite('Pnewfile.txt',P_out);
    dlmwrite('Nnewfile.txt',N_out);
    dlmwrite('Enewfile.txt',E_out);
    dlmwrite('Mnewfile.txt',M_out);    
    
    concat(Pfilename,'Pnewfile.txt')
    concat(Nfilename,'Nnewfile.txt')
    concat(Efilename,'Enewfile.txt')
    concat(Mfilename,'Mnewfile.txt')
    
    % set new initial conditions
    y0 = y_out(size(ts_vec,1),:);
        
end



% plot initial & final distributions
%     figure
%     plot((1:1:Pdim1),P0)
%     
%     figure
%     Pfin = squeeze(P_out(n_ts,:));
%     plot((1:1:Pdim1),Pfin)
    
    figure
    hold on
    surf(P_out,'MeshStyle','row')
    hold off
    axis([0 Pdim1 0 n_ts 0 max(P_out(:,x0))])
%    set(gca,'ZScale','log')

    figure
    hold on
    surf(N_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 N0density])

    figure
    hold on
    surf(E_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 5000])
    
    figure
    hold on
    surf(M_out,'MeshStyle','row')
    hold off
    axis([0 Ldim1 0 n_ts 0 400])



    Ptot = sum(P_out,2);
    Ntot = sum(N_out,2);
    Etot = sum(E_out,2);
    Mtot = sum(M_out,2);
    figure
    semilogy(ts_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
%    axis([0 days 1 10^8])
