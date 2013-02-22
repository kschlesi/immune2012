function dy = ss_dy(t,y,Pdim1,Ldim1)

global r_ h_ sigma_ de_ f_ k_ c gammas1D lambdas1D mu_ R_ dh_ ;

% create separate P, N, E, M vectors
P = y(1:Pdim1);
N = y(Pdim1+1:Pdim1+Ldim1);
E = y(Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M = y(Pdim1+2*Ldim1+1:end);

% create stochastic mutation matrix (size Pdim1 x Pdim1)
%mrates = eye(Pdim1);
mrates = zeros(Pdim1,Pdim1);
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/((i-j)^c));
            mrates(j,i) = mrates(i,j);
        end
        iloss = sum(mrates(i,:));
        mrates(i,i) = 1-iloss;
    end
        iloss = sum(mrates(1,:));
        mrates(1,1) = 1-iloss;
    
    
% enforcing P cutoff for calculating everything...
for i=1:Pdim1
    if(P(i)<mu_)
        P(i)=0;
    end
    if(N(i)<mu_)
        N(i)=0;
    end
%     if(E(i)<mu_)
%         E(i)=0;
%     end
%     if(M(i)<mu_)
%         M(i)=0;
%     end
end
        
% calculate dP (all size Pdim1 x 1)        
dmut = zeros(Pdim1,1);
omega = zeros(Pdim1,1);
    for i=1:Pdim1
        dmut(i) = r_.*squeeze(sum(P.*squeeze(mrates(:,i))));
        omega(i) = sum(shiftdim(gammas1D(i,:)).*(N + M + E));
    end
dP = dmut.*lambdas1D - h_.*omega.*P;
ndP = 0;
for i=1:Pdim1
    if(dP(i)<mu_)
        dP(i)=0;
        ndP = ndP+1000
    end
end

% calculate dL's (all size Ldim1 x 1)
Pofy = zeros(Ldim1,1);
    for j = 1:Ldim1
        Pofy(j)= sum(P.*squeeze(gammas1D(:,j)));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy);
dN = -sigma_.*N.*satfunc;
dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,1)-satfunc);
dM = f_.*de_.*E.*(ones(Ldim1,1)-satfunc) - sigma_.*M.*satfunc;

t+ndP

dy = [dP;dN;dE;dM];
