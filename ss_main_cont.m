% fitness landscape_cont

clear

global r_ h_ sigma_ de_ f_ k_ c b beta_ p_ ;

days = 0;
stepsize = 0.1;
olddays = 5;
oldss = 0.1;

% file to which new days will be appended
Pfilename = 'Ploop5.txt';
Nfilename = 'Nloop5.txt';
Efilename = 'Eloop5.txt';
Mfilename = 'Mloop5.txt';

% dimensions of 1D shape space
Pdim1 = 600;
Ldim1 = 600;
x0 = 300;
p_ = (1-exp(-1*((Pdim1)^2)/(8*beta_^2)));

% gammas & lambdas
gammas1D = zeros(Pdim1,Ldim1);
lambdas1D = zeros(Pdim1,1);
for i=1:Pdim1;
    lambdas1D(i) = p_*(1-exp(-1*((i-x0)^2)/(2*beta_^2)));
    for j=1:Ldim1;
        gammas1D(i,j) = exp(-1*((i-j)^2)/(2*b^2));
    end
end

% initial inoculation in shape space
P0in = csvread(Pfilename);
N0in = csvread(Nfilename);
E0in = csvread(Efilename);
M0in = csvread(Mfilename);
old_ts = size(P0in,1);
P0 = shiftdim(P0in(old_ts,:),1);
N0 = shiftdim(N0in(old_ts,:),1);
E0 = shiftdim(E0in(old_ts,:),1);
M0 = shiftdim(M0in(old_ts,:),1);


% creating initial conditions vector
y0 = [P0;N0;E0;M0];

% plotting initial conitions
    % figure
    % semilogy((1:Pdim1+3*Ldim1),(shiftdim(y0)))

% integrating diffeqs in time with a FOR LOOP
options = odeset('AbsTol',1e-3);
nsteps = cast(days/stepsize,'uint16');
n_ts = old_ts;
n_ts
for j=1:nsteps
    
    % integrate between two external steps (of size stepsize)
    i = cast(j,'double');
    [ts_vec,y_out] = ode45(@(t,y)ss_dy(t,y,Pdim1,Ldim1),[(i-1)*stepsize,i*stepsize],y0,options);

    % add new internal steps to overall n_ts
    size(ts_vec,1)
    n_ts = n_ts + size(ts_vec,1)-1;

    % implement one-cell cutoff for all P
    for pcount=1:Pdim1
        if(y_out(end,pcount)<1)
            y_out(end,pcount)=0;
        end
    end
        
   % save & append y-output
    P_out = y_out(:,1:Pdim1);
    N_out = y_out(:,Pdim1+1:Pdim1+Ldim1);
    E_out = y_out(:,Pdim1+Ldim1+1:Pdim1+2*Ldim1);
    M_out = y_out(:,Pdim1+2*Ldim1+1:end);

    dlmwrite('Pnewfile.txt',P_out);
    dlmwrite('Nnewfile.txt',N_out);
    dlmwrite('Enewfile.txt',E_out);
    dlmwrite('Mnewfile.txt',M_out);    
    
    concat(Pfilename,'Pnewfile.txt')
    concat(Nfilename,'Nnewfile.txt')
    concat(Efilename,'Enewfile.txt')
    concat(Mfilename,'Mnewfile.txt')
    
    % set new initial conditions
    y0 = y_out(size(ts_vec,1),:);
        
end


% plot initial & final distributions
    figure
    plot((1:1:Pdim1),P0)
    
    figure
    Pfin = squeeze(P_out(end,:));
    plot((1:1:Pdim1),Pfin)
    
    Ptot = sum(P_out,2);
    Ntot = sum(N_out,2);
    Etot = sum(E_out,2);
    Mtot = sum(M_out,2);
    figure
    semilogy(ts_vec,Ptot,ts_vec,Ntot+Mtot+Etot)
