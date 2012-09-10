% ShSpSimdy function handle for ode45 integration routine
function dy = ShSpSimdy(t,y,dim1,dim2)

global b_ r_ h_ sigma_ k_ de_ f_ gammas_ ;

% separate cell types, leaving everything linear
% each NEM variable vector has length dim1*dim2
P = y(1);
Nlin = y(2:dim1*dim2+1);
Elin = y(dim1*dim2+2:2*dim1*dim2+1);
Mlin = y(2*dim1*dim2+2:end);
gammas_lin = reshape(gammas_,dim1*dim2,1); %same vec lgth as NEMs

% calculating dCells with diffeqs (linear cell vectors)
    satfunc = (P/(k_+P));
    omega_ = sum(gammas_lin.*(Nlin+Mlin+Elin)); 
    dP = r_ * P - h_ * P* omega_;
    dN = -sigma_ * gammas_lin .* Nlin *satfunc;
    dE = 2 * sigma_ * gammas_lin .* Nlin * satfunc + sigma_ * gammas_lin .* Elin * satfunc - de_ * Elin * (1 - satfunc);
    dM = f_ .* de_ .* Elin .* (1-satfunc);
    
dy = [dP;dN;dE;dM];



