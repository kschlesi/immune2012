function dy = ss_dy(t,y,b0,gammas1D,lambdas1D)

% parse b0 and set constant parameters in scope
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
%chi_ = b0(16);
Qstep = b0(17);
Gamma_ = b0(18);
%nrandon = b0(19);
delta_ = b0(20);
%muton = b0(21);
spliton = b0(22);
%pinit = b0(23);
%nu_ = b0(24);

% create separate P, L vectors
% & set all < mu_ pops to 0 in ss_dy only (not returned to ss_main)
P = y(1:Pdim1);
L = y(Pdim1+1:Pdim1+Ldim1);
P = P.*(P>=mu_);
L = L.*(L>=mu_);
Pis0 = (P<mu_);  % keep track of whether each site is below mu_

% calculate dP (all size Pdim1 x 1)        
% PRmat = repmat(P.*lambdas1D,1,Pdim1);
% dmut = sum(PRmat.*mrates,1)';
% clear PRmat;
Lmat = repmat(L',Ldim1,1);
omega = sum(Lmat.*gammas1D,2);
clear Lmat;
Ptot = sum(P);
dP = (r_.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega).*P;
Pmat = repmat(P,1,Pdim1);
if spliton
    %dmut = sum(Pmat.*mrates,1)';
    dP = (r_.*lambdas1D.*(1-capon*Ptot/K_) - h_.*omega).*P;% + dmut;
end
zerodP = Pis0.*(dP<(mu_/Qstep)); % zero dP if Pis0, unless dP > mu_ (per mutation step)
%ndP = sum(zerodP);               % number of sites that were at 0, and stayed there
dP = dP.*(1-zerodP); 

% calculate dL's (all size Ldim1 x 1)
if hsaton
    Hsat = (sum(L) - R_); % for lymphocyte constraint
else
    Hsat = 0; % if no constraint
end
Pofy = sum(Pmat.*gammas1D)';
clear Pmat;
satfunc = Pofy./(k_.*ones(Ldim1,1)+Pofy); % pathogen saturation function
dL = Gamma_ + (sigma_*satfunc - delta_*(1-satfunc) - dh_*Hsat).*L;

% tstep = t-tlast;
% if muton*tstep>0
%     ploss = poissrnd(P.*nu_.*tstep);
%     newsites = floor(exprnd(chi_,[size(ploss,1),max(ploss)])+1);
%     for i=1:size(ploss,1)
%         newsites(i,ploss(i)+1:end) = zeros;
%     end
%     mrates = histc(newsites',linspace(1,Pdim1,Pdim1))';
%     pgain = sum(mrates,1)';
% end

%tlast = t;

% print to command window: time (in days) and 
% # of sites with NO pathogen both before and after this timestep
disp(t);
%disp(ndP);

dy = [dP;dL];

end