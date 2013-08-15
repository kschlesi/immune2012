function lambas = Lambdas(eps_,Pdim1)

disp(eps_);

if ~eps_
    lambas = ones(Pdim1,1);
    lambas(1:2) = [0;0.5];
    lambas(Pdim1-1:end) = [0.5;0];
else
    lambas = 1 - (2*eps_)./(Pdim1 + 2*eps_ - abs(Pdim1-2*(1:Pdim1)'));
end

end