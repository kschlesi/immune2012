function dy = ss_dy(t,y,Pdim1,Ldim1)

global r_ h_ sigma_ de_ f_ k_ c gammas1D lambdas1D mu_ R_ dh_ K_ Gamma_ ;
global chi_ Qstep tgone Nstep ntgone nrandon capon hsaton mrates ;

% create separate P, N, E, M vectors
P = y(1:Pdim1);
N = y(Pdim1+1:Pdim1+Ldim1);
E = y(Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M = y(Pdim1+2*Ldim1+1:end);

% create stochastic mutation matrix (size Pdim1 x Pdim1)
% (matrix is symmetric, and all rows sum to 1)
% mrates = eye(Pdim1);      % use this line for no mutation
if (t-tgone)>=Qstep         % use this 'if statement' for mutation every Qstep
    mrates = zeros(Pdim1,Pdim1);
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/(i-j)^c)/chi_; % mutation probability
            mrates(j,i) = mrates(i,j);
        end
        iloss = sum(mrates(i,:));
        mrates(i,i) = 1-iloss;
    end
    iloss = sum(mrates(1,:));
    mrates(1,1) = 1-iloss;
    tgone = t;   
end
    
% enforcing cutoff for calculating everything...
% this sets all < mu_ pops to 0 in ss_dy only (not returned to ss_main)
Pis0 = zeros(Pdim1,1);  % keep track of whether each site is below mu_
Nis0 = zeros(Ldim1,1);
Eis0 = zeros(Ldim1,1);
Mis0 = zeros(Ldim1,1);
for i=1:Pdim1
    if(P(i)<mu_)
        P(i)=0;
        Pis0(i)=1;
    end
end
for i=1:Ldim1
    if(N(i)<mu_)
        N(i)=0;
        Nis0(i) = 1;
    end
    if(E(i)<mu_)
        E(i)=0;
        Eis0(i) = 1;
    end
    if(M(i)<mu_)
        M(i)=0;
        Mis0(i) = 1;
    end
end
        
% calculate dP (all size Pdim1 x 1)        
dmut = zeros(Pdim1,1);
omega = zeros(Pdim1,1);
    for i=1:Pdim1
        dmut(i) = squeeze(sum(P.*squeeze(mrates(:,i)))); % mut matrix mrates*P & summed
        omega(i) = sum(shiftdim(gammas1D(i,:)).*(N + M + E)); % effectivity
    end
Ptot = sum(P);
dP = r_.*dmut.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega.*P;
ndP = 0;
for i=1:Pdim1         %% IF Pis0 (that is, we COUNT no P there, or P < mu_)  
    if(Pis0(i)==1 && dP(i)<(mu_/Qstep))  %% THEN P cannot show up there (dP = 0)
        dP(i)=0;                         %% UNLESS dP > mu_ (per mutation step)
        ndP = ndP+1; % ndP = number of sites that were at 0,
    end              % and received mutations too small
end                  % to overcome cutoff

% calculate dL's (all size Ldim1 x 1)
if hsaton
    Hsat = (sum(N + E + M) - R_); % for lymphocyte constraint
else
    Hsat = 0; % if no constraint
end
Pofy = zeros(Ldim1,1); % total P weighted by affinity to lymphocyte
    for j = 1:Ldim1
        Pofy(j)= sum(P.*squeeze(gammas1D(:,j)));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy); % pathogen saturation function
Nflux = Gamma_; % naive cell influx
if(nrandon)
    if (ntgone-t)>=Nstep             %% re-distributes incoming naive cells
        Nflux = unifrndpop(Ldim1,Gamma_,mu_);    %% at intervals of 'Nstep'
        ntgone = t;
    end
end
dN = Nflux - sigma_.*N.*satfunc - dh_.*Hsat.*N;
dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,1)-satfunc) - dh_.*Hsat.*E;
dM = f_.*de_.*E.*(ones(Ldim1,1)-satfunc) - sigma_.*M.*satfunc - dh_.*Hsat.*M;

% prints to command window: time (in days) and # of sites with NO pathogen
t
ndP

dy = [dP;dN;dE;dM];