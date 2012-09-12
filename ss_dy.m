function dy = ss_dy(t,y,Pdim1,Ldim1)

global r_ h_ sigma_ de_ f_ k_ c gammas1D lambdas1D;

% create separate P, N, E, M vectors
P = y(1:Pdim1);
N = y(Pdim1+1:Pdim1+Ldim1);
E = y(Pdim1+Ldim1+1:Pdim1+2*Ldim1);
M = y(Pdim1+2*Ldim1+1:end);

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
	    % enforcing one-cell minimum at each pathogen location
	    if (P(i,j)<1)
	        P(i,j)=0;
	    end
	end
    end

    dP = r_.*P = h_.*P.*omega_;

    % all below are size = Ldim1,Ldim2
    for i = 1:Ldim1
        for j = 1:Ldim2
            Pofy(i,j)= sum(sum(P.*squeeze(gamma_lib(:,:,i,j)),1),2);
        end
    end

satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy);
dN = -sigma_.*N.*satfunc;
dE = sigma_.*(2*N + E + 2*M).*satfunc - de_.*E.*(ones(Ldim1,1)-satfunc);
dM = f_.*de_.*E.*(ones(Ldim1,1)-satfunc) - sigma_.*M.*satfunc;

t

dy = [dP;dN;dE;dM];
