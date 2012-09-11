% P_eff function handle for ode45 integration routine
function dy = ss_dy(t,y,Pdim1,Pdim2,Ldim1,Ldim2)

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
Pofy = zeros(Ldim1,Ldim2);

% calculating dCells with diffeqs
    for i = 1:Pdim1
        for j = 1:Pdim2
            omega_(i,j) = sum(sum(squeeze(gamma_lib(i,j,:,:)).*( N + M + E ),1),2); % size = Pdim1,Pdim2
        end
    end
    dP = r_.*P - h_.*P.*omega_;
    % all below are size = Ldim1,Ldim2
    for i = 1:Ldim1
        for j = 1:Ldim2
            Pofy(i,j)= sum(sum(P.*squeeze(gamma_lib(:,:,i,j)),1),2);
        end
    end
    satfunc = Pofy./(k_.*ones(Ldim1,Ldim2)+Pofy);
    dN = -sigma_.*N.*satfunc;
    dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,Ldim2)-satfunc);
    dM = f_.*de_.*E.*(ones(Ldim1,Ldim2)-satfunc)-sigma_.*M.*satfunc;

t
N(30,25)
satfunc(30,25)
    
% reshaping and creating final dy vector
dPlin = reshape(dP,Pdim1*Pdim2,1);
dNlin = reshape(dN,Ldim1*Ldim2,1);
dElin = reshape(dE,Ldim1*Ldim2,1);
dMlin = reshape(dM,Ldim1*Ldim2,1);
dy = [dPlin;dNlin;dElin;dMlin];
function dy = Pmutn(t,y,Pdim1,Pdim2)

global r_ mrate ijs;

P = reshape(y,Pdim1,Pdim2);

% random walk steps of each location, added to location itself, for matrix
% of new locations
rwsteps = cast(unifrnd(-1.5,1.5,Pdim1,Pdim2,2),'int16') + ijs;

Pmults = zeros(Pdim1,Pdim2); % scalar field of #pathogen mutating INTO each P site
for i=0:Pdim1-1
    for j=0:Pdim2-1
        if squeeze(rwsteps(mod(i+1,Pdim1)+1,j+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(mod(i+1,Pdim1)+1,j+1);
        end
        if squeeze(rwsteps(mod(i+Pdim1-1,Pdim1)+1,j+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(mod(i+Pdim1-1,Pdim1)+1,j+1);
        end
        if squeeze(rwsteps(i+1,mod(j+1,Pdim2)+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(i+1,mod(j+1,Pdim2)+1);
        end
        if squeeze(rwsteps(i+1,mod(j+Pdim2-1,Pdim2)+1,:))==[i+1;j+1]
            Pmults(i+1,j+1) = Pmults(i+1,j+1)+P(i+1,mod(j+Pdim2-1,Pdim2)+1);
        end
    end
end

t

dP = mrate.*(Pmults-P)+r_.*P;
dy = reshape(dP,Pdim1*Pdim2,1);
