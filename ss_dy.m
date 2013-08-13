function dy = ss_dy(t,y,b0,gammas1D,lambdas1D)

global mrates tgone ;

% varrs = {'r_';'h_';'sigma_';'k_';'b';'eps_';'mu_';...
%     'dh_';'K_';'R_';'capon';'hsaton';'Pdim1';'Ldim1';'x0';'chi_';'Qstep';...
%     'Gamma_';'nrandon';'delta_';'muton';'pinit'};
% for i=1:size(b0,1)
%     eval([char(varrs{i,1}) ' = ' num2str(b0(i,1)) ';']);
% end
r_ = b0(1);
h_ = b0(2);
sigma_ = b0(3);
k_ = b0(4);
%b = b0(5);
%eps_ = b0(6);
mu_ = b0(7);
dh_ = b0(8);
K_ = b0(9);
R_ = b0(10);
capon = b0(11);
hsaton = b0(12);
Pdim1 = b0(13);
Ldim1 = b0(14);
%x0 = b0(15);
chi_ = b0(16);
Qstep = b0(17);
Gamma_ = b0(18);
%nrandon = b0(19);
delta_ = b0(20);
muton = b0(21);
%pinit = b0(22);


% create separate P, L vectors
P = y(1:Pdim1);
L = y(Pdim1+1:Pdim1+Ldim1);

% create stochastic mutation matrix (size Pdim1 x Pdim1)
% (matrix is symmetric, and all rows sum to 1)
if (muton)              % use this 'if statement' for mutation every Qstep
    if (t-tgone)>=Qstep         
        mrates = Qmatrix(Pdim1,chi_);
        tgone = t;   
    end
end    

% enforcing cutoff for calculating everything...
% this sets all < mu_ pops to 0 in ss_dy only (not returned to ss_main)
P = P.*(P>=mu_);
L = L.*(L>=mu_);
Pis0 = ones(Pdim1,1)-(P>=mu_);  % keep track of whether each site is below mu_
        
% calculate dP (all size Pdim1 x 1)        
dmut = zeros(Pdim1,1);
omega = zeros(Pdim1,1);
    for i=1:Pdim1
        dmut(i) = sum(P.*mrates(:,i)); % mut matrix mrates*P & summed
        omega(i) = sum(gammas1D(:,i).*L); % effectivity
    end
Ptot = sum(P);
dP = (r_.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega).*P + dmut;

% IF Pis0 (that is, P < mu_) THEN P cannot show up there (dP = 0)
% UNLESS dP > mu_ (per mutation step);
% ndP = number of sites that were at 0, and stayed there
zeroP = Pis0.*(dP<(mu_/Qstep));
ndP = sum(zeroP);
dP = dP.*(1-zeroP);

% calculate dL's (all size Ldim1 x 1)
if hsaton
    Hsat = (sum(L) - R_); % for lymphocyte constraint
else
    Hsat = 0; % if no constraint
end
Pofy = zeros(Ldim1,1); % total P weighted by affinity to lymphocyte
    for j = 1:Ldim1
        Pofy(j)= sum(P.*gammas1D(:,j));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy); % pathogen saturation function
dL = Gamma_.*ones(Ldim1,1) + (sigma_.*satfunc - delta_.*(ones(Ldim1,1)-satfunc) - dh_.*Hsat.*ones(Ldim1,1)).*L;

% print to command window: time (in days) and # of sites with NO pathogen
disp(t);
disp(ndP);

dy = [dP;dL];

end