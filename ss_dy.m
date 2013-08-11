function dy = ss_dy(t,y,b0,gammas1D,lambdas1D)

global mrates tgone ;

varrs = {'r_';'h_';'sigma_';'k_';'b';'eps_';'mu_';...
    'dh_';'K_';'R_';'capon';'hsaton';'Pdim1';'Ldim1';'x0';'chi_';'Qstep';...
    'Gamma_';'nrandon';'delta_';'muton';'pinit'};
for i=1:size(b0,1)
    eval([char(varrs{i,1}) ' = ' num2str(b0(i,1))]);
end
%tgone = tcounts(1);
%ntgone = tcounts(2);

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
Pis0 = zeros(Pdim1,1);  % keep track of whether each site is below mu_
Lis0 = zeros(Ldim1,1);
for i=1:Pdim1
    if(P(i)<mu_)
        P(i)=0;
        Pis0(i)=1;
    end
end
for i=1:Ldim1
    if(L(i)<mu_)
        L(i)=0;
        Lis0(i) = 1;
    end
end
        
% calculate dP (all size Pdim1 x 1)        
dmut = zeros(Pdim1,1);
omega = zeros(Pdim1,1);
    for i=1:Pdim1
        dmut(i) = squeeze(sum(P.*squeeze(mrates(:,i)))); % mut matrix mrates*P & summed
        omega(i) = sum(shiftdim(gammas1D(i,:)).*L); % effectivity
    end
Ptot = sum(P);
dP = r_.*dmut.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega.*P;
ndP = 0;
for i=1:Pdim1         %% IF Pis0 (that is, we COUNT no P there, or P < mu_)  
    %if(Pis0(i)==1 && dP(i)<(mu_/Qstep))  %% THEN P cannot show up there (dP = 0)
    if(Pis0(i)==1 && dP(i)<(mu_/0.1))  
        dP(i)=0;                         %% UNLESS dP > mu_ (per mutation step)
        ndP = ndP+1; % ndP = number of sites that were at 0,
    end              % and received mutations too small
end                  % to overcome cutoff

% calculate dL's (all size Ldim1 x 1)
if hsaton
    Hsat = (sum(L) - R_); % for lymphocyte constraint
else
    Hsat = 0; % if no constraint
end
Pofy = zeros(Ldim1,1); % total P weighted by affinity to lymphocyte
    for j = 1:Ldim1
        Pofy(j)= sum(P.*squeeze(gammas1D(:,j)));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy); % pathogen saturation function
dL = Gamma_.*ones(Ldim1,1) + (sigma_.*satfunc - delta_.*(ones(Ldim1,1)-satfunc) - dh_.*Hsat.*ones(Ldim1,1)).*L;

% prints to command window: time (in days) and # of sites with NO pathogen
disp(t);
disp(ndP);

dy = [dP;dL];

end