function lambdas = Lambdas(eps_,Pdim1)

if ~eps_
    lambdas = ones(Pdim1,1);
    lambdas(1:2) = [0;0.5];
    lambdas(Pdim1-1:end) = [0.5;0];
else
    lambdas = zeros(Pdim1,1);           
    for i=1:Pdim1
        lambdas = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*i));
    end
end

end