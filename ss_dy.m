function dy = ss_dy(t,y,Pdim1,Ldim1)

global r_ h_ sigma_ de_ f_ k_ c gammas1D lambdas1D mu_ R_ dh_ K_ ;
global chi_ Qstep tgone Nstep ntgone nrandon capon hsaton mrates ;

% create separate P, N, E, M vectors
P = y(1:Pdim1);
N = y(Pdim1+1:Pdim1+Ldim1);
E = y(Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M = y(Pdim1+2*Ldim1+1:end);

% create stochastic mutation matrix (size Pdim1 x Pdim1)
% mrates = eye(Pdim1);      % use this line for no mutation
if (t-tgone)>=Qstep         % use this if statement for mutation every Qstep
    mrates = zeros(Pdim1,Pdim1);
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/(i-j)^c)/chi_;
            mrates(j,i) = mrates(i,j);
        end
        iloss = sum(mrates(i,:));
        mrates(i,i) = 1-iloss;
    end
    iloss = sum(mrates(1,:));
    mrates(1,1) = 1-iloss;
    tgone = t;   
end
    
% enforcing P cutoff for calculating everything...
Pis0 = zeros(Pdim1,1);
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
        dmut(i) = squeeze(sum(P.*squeeze(mrates(:,i))));
        omega(i) = sum(shiftdim(gammas1D(i,:)).*(N + M + E));
    end
Ptot = sum(P);
dP = r_.*dmut.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega.*P;
ndP = 0;
for i=1:Pdim1   %% IF Pis0 (that is, we COUNT no P there, or P < mu_)  
    if(Pis0(i)==1 && dP(i)<(mu_/Qstep))  %% THEN P cannot show up there (dP = 0)
        dP(i)=0;                 %% UNLESS dP > mu_ (permutation... may want to change)
        ndP = ndP+1;
    end
end

% calculate dL's (all size Ldim1 x 1)
if hsaton
    Hsat = (sum(N + E + M) - R_);
else
    Hsat = 0;
end
Pofy = zeros(Ldim1,1);
    for j = 1:Ldim1
        Pofy(j)= sum(P.*squeeze(gammas1D(:,j)));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy);
Nflux = Gamma_;
if(nrandon)
    if (ntgone-t)>=Nstep
        Nflux = unifrndpop(Ldim1,Gamma_,mu_);
        ntgone = t;
    end
end
dN = Nflux - sigma_.*N.*satfunc - dh_.*Hsat.*N;
dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,1)-satfunc) - dh_.*Hsat.*E;
dM = f_.*de_.*E.*(ones(Ldim1,1)-satfunc) - sigma_.*M.*satfunc - dh_.*Hsat.*M;

% ndN = 0;
% ndE = 0;
% ndM = 0;
% for i=1:Ldim1   %% IF Nis0 (that is, we COUNT no N there, or N < mu_)  
%     if(Nis0(i)==1 && dN(i)<(mu_/Qstep))  %% THEN N cannot show up there (dN = 0)
%         dN(i)=0;                 %% UNLESS dN > mu_ (permutation... may want to change)
%         ndN = ndN+1;
%     end
%                  %% IF Eis0 (that is, we COUNT no E there, or E < mu_)  
%     if(Eis0(i)==1 && dE(i)<(mu_/Qstep))  %% THEN E cannot show up there (dE = 0)
%         dE(i)=0;                 %% UNLESS dE > mu_ (permutation... may want to change)
%         ndE = ndE+1;
%     end
%                  %% IF Mis0 (that is, we COUNT no M there, or M < mu_)  
%     if(Mis0(i)==1 && dM(i)<(mu_/Qstep))  %% THEN M cannot show up there (dM = 0)
%         dM(i)=0;                 %% UNLESS dM > mu_ (permutation... may want to change)
%         ndM = ndM+1;
%     end
% end




t
ndP

dy = [dP;dN;dE;dM];
