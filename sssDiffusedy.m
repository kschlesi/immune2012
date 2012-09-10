% ShSpSimdy function handle for ode45 integration routine
function dy = sssDiffusedy(t,y,Pdim1,Pdim2,Ldim1,Ldim2)

global r_ h_ sigma_ k_ de_ f_ gamma_lib ;

% separate cell types
Plin = y(1:Pdim1*Pdim2);
Nlin = y(Pdim1*Pdim2+1:Pdim1*Pdim2+Ldim1*Ldim2);
Elin = y(Pdim1*Pdim2+Ldim1*Ldim2+1:Pdim1*Pdim2+2*Ldim1*Ldim2);
Mlin = y(Pdim1*Pdim2+2*Ldim1*Ldim2+1:end);

P = reshape(Plin,Pdim1,Pdim2);
N = reshape(Nlin,Ldim1,Ldim2);
E = reshape(Elin,Ldim1,Ldim2);
M = reshape(Mlin,Ldim1,Ldim2);

% hold on
% figure
% surf(P)
% hold off

% preallocating arrays
omega_ = zeros(Pdim1,Pdim2);
dN = zeros(Ldim1,Ldim2);
dE = zeros(Ldim1,Ldim2);
dM = zeros(Ldim1,Ldim2);

% calculating dCells with diffeqs
    satfunc = (1/(Pdim1*Pdim2)).*P./(k_.*ones(Pdim1,Pdim2)+P); % size = Pdim1,Pdim2, normalised for summing already
    for i = 1:Pdim1
        for j = 1:Pdim2
            omega_(i,j) = sum(sum(squeeze(gamma_lib(i,j,:,:)).*( N + M + E ),1),2); % size = Pdim1,Pdim2
        end
    end
    dP = r_.*P - h_.*P.*omega_;
    % all below are size = Ldim1,Ldim2
    for i = 1:Ldim1
        for j = 1:Ldim2
            dN(i,j)= -sigma_*N(i,j)*sum(sum(satfunc.*squeeze(gamma_lib(:,:,i,j)),1),2);
            dE(i,j)= sigma_*(2*N(i,j)+E(i,j)+2*M(i,j))*sum(sum(satfunc.*squeeze(gamma_lib(:,:,i,j)),1),2) - de_*E(i,j)*(1-sum(sum(satfunc,1),2));
            dM(i,j)= f_*de_*E(i,j)*(1-sum(sum(satfunc,1),2))-sigma_*M(i,j)*sum(sum(satfunc.*squeeze(gamma_lib(:,:,i,j)),1),2);
        end
    end

t
r_
h_*omega_(30,25)
N(30,25)
sigma_*sum(sum(satfunc.*squeeze(gamma_lib(:,:,30,25)),1),2)
    
% reshaping and creating final dy vector
dPlin = reshape(dP,Pdim1*Pdim2,1);
dNlin = reshape(dN,Ldim1*Ldim2,1);
dElin = reshape(dE,Ldim1*Ldim2,1);
dMlin = reshape(dM,Ldim1*Ldim2,1);

dy = [dPlin;dNlin;dElin;dMlin];
