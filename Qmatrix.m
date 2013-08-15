function mrates = Qmatrix(Pdim1,chi_,spliton)

notnormz = 1;
iminj = repmat((1:Pdim1)',1,Pdim1)-repmat(1:Pdim1,Pdim1,1)+eye(Pdim1);

if ~chi_
    mrates = zeros(Pdim1);
    if ~spliton
        % uncomment below if combined mutation and growth terms
        mrates = eye(Pdim1);
    end
else
    while notnormz
        notnormz = 0;
        mrates = (1/(Pdim1*chi_)).*abs(randn(Pdim1))./(iminj.^2).*(triu(ones(Pdim1),1));
        mrates = mrates + triu(mrates,1)';
        iloss = sum(mrates,1);   % if unsymmetrising, check this!!!!
        notnormz = sum((iloss>1));
    end
    if ~spliton
        % uncomment below if combined mutation and growth terms
        mrates = mrates + diag(1-iloss);
    end
end
    

end