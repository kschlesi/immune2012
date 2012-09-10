function dy = ss_dy(t,y,Pdim1,Ldim1)

global r_ h_ sigma_ de_ f_ k_ D gammas1D;

% create separate P, N, E, M vectors
P = y(1:Pdim1);
N = y(Pdim1+1:Pdim1+Ldim1);
E = y(Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M = y(Pdim1+2*Ldim1+1:end);


dmut = zeros(Pdim1,1);
omega = zeros(Pdim1,1);
    for i=2:Pdim1-1
        dmut(i) = D*(P(i-1)-2*P(i)+P(i+1));
        omega(i) = sum(shiftdim(gammas1D(i,:)).*(N + M + E));
    end
    % on boundaries, use forward/backward differences
    dmut(1) = D*(P(1)-2*P(2)+P(3));
    dmut(Pdim1) = D*(P(Pdim1)-2*P(Pdim1-1)+P(Pdim1-2));
    omega(1) = sum(shiftdim(gammas1D(1,:)).*(N + M + E));
    omega(Pdim1) = sum(shiftdim(gammas1D(Pdim1,:)).*(N + M + E));
    
%dP = dmut;
dP = dmut + r_.*P - h_.*omega.*P;

Pofy = zeros(Ldim1,1);
    for j = 1:Ldim1
        Pofy(j)= sum(P.*squeeze(gammas1D(:,j)));
    end
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy);
dN = -sigma_.*N.*satfunc;
dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,1)-satfunc);
dM = f_.*de_.*E.*(ones(Ldim1,1)-satfunc)-sigma_.*M.*satfunc;

t

dy = [dP;dN;dE;dM];
