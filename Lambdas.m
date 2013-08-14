function lambdas = Lambdas(eps_,Pdim1)

disp(eps_);

if ~eps_
    lambdas = ones(Pdim1,1);
    lambdas(1:2) = [0;0.5];
    lambdas(Pdim1-1:end) = [0.5;0];
else
    disp(size((Pdim1 + 2*eps_ - abs(Pdim1-2*(1:Pdim1)'))));
    lambdas = 1 - (2*eps_)/(Pdim1 + 2*eps_ - abs(Pdim1-2*(1:Pdim1)'));
    figure
    plot(lambdas)
end

end