% shape space simulation 1
% with loop for random walk in shapespace

clear

global b_ r_ h_ sigma_ k_ de_ f_ gammas_ j ;

N0density = 3; 
P0 = 10; 
r_ = 3.3;
h_ = 10^-5;
sigma_ = 2;
k_ = 10^5;
de_ = 0.35;
f_ = 0.1;
b_ = 5;
days = 20;

% dimensions of 2D shape space
dim1 = 60;
dim2 = 50;

% initial location of inoculation in shape space
x0 = [30,25];

% creating initial gammas and lymphocytes
N0 = N0density.*ones(dim1,dim2);
E0 = zeros(dim1,dim2);
M0 = zeros(dim1,dim2);
N0lin = reshape(N0,dim1*dim2,1);
E0lin = reshape(E0,dim1*dim2,1);
M0lin = reshape(M0,dim1*dim2,1);
gammas_ = Gammas(x0,N0,1,b_);

% ::mutation-step loop for randomwalking antigen::

% setting loop parameters:
m_perday = 10;          % FIX TYPES PROPERLY
day_step = 1/m_perday;
stepsize = 0.01;
n_msteps = cast(days*m_perday,'uint16');  % should be int
n_ts = cast(day_step/stepsize,'uint16');  % should be int

% preallocating memory for storage arrays
xstep = zeros(n_msteps,2);     % stores antigen location at mstep
xstep(1,:) = x0;
y_out_mstep = zeros(n_msteps,n_ts,3*dim1*dim2+1); % stores ode45 output at mstep


% mutating infection loop
for i = 1:n_msteps
    
    if (mod(i,15)==0)
    % having set xstep(i) & calculated gammas... 
        % plotting new gammas
        hold on
        figure
        surf(gammas_)
        hold off
    end;

    % creating vec of initial conditions
%   N0lin = N0lin + squeeze(M_outlin(n_ts,:));
    y0 = [P0;N0lin;E0lin;M0lin];

    % integrating diffeqs
    j = cast(i,'double');
    [ts_vec,y_out] = ode45(@(t,y)ShSpSimdy(t,y,dim1,dim2),(day_step*(j-1):stepsize:(day_step*j)-stepsize),y0);
%   n_ts = size(ts_vec,1);
    y_out_mstep(i,:,:) = y_out;
        
    % set new location for virus (random walk in shapespace)
    newxx = xstep(i,1)+randn(1);
    newxy = xstep(i,2)+randn(1);
    xstep(i+1,:) = [newxx,newxy];
    
    % creating gammas for new location
    gammas_ = Gammas(xstep(i+1,:),N0,1,b_);

    % overwrite N0, E0, M0
    P0 = y_out(n_ts,1);
    N0lin = y_out(n_ts,2:dim1*dim2+1)';
    E0lin = y_out(n_ts,dim1*dim2+2:2*dim1*dim2+1)';
    M0lin = y_out(n_ts,2*dim1*dim2+2:end)';

end

    % setting final values
    P_o = y_out_mstep(:,:,1);
    N_o = y_out_mstep(:,:,2:dim1*dim2+1);
    E_o = y_out_mstep(:,:,dim1*dim2+2:2*dim1*dim2+1);
    M_o = y_out_mstep(:,:,2*dim1*dim2+2:end);
    
    P_out = reshape(P_o',n_msteps*n_ts,1);
    N_out = reshape(N_o,n_msteps,n_ts,dim1,dim2);
    E_out = reshape(E_o,n_msteps,n_ts,dim1,dim2);
    M_out = reshape(M_o,n_msteps,n_ts,dim1,dim2);

    N_tot = reshape(sum(sum(N_out,4),3)',n_msteps*n_ts,1); 
    E_tot = reshape(sum(sum(E_out,4),3)',n_msteps*n_ts,1);
    M_tot = reshape(sum(sum(M_out,4),3)',n_msteps*n_ts,1);

    t_vec = (0:stepsize:size(P_out)*stepsize-stepsize);    %day_step*j-stepsize]';
    
    % plotting immune response
    figure
    semilogy(t_vec,P_out,t_vec,(M_tot + N_tot + E_tot))%,t_vec,M_tot,t_vec,N_tot,t_vec,E_tot)
    axis([0 days 1 10^12])

    
