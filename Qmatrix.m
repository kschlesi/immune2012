function mrates = Qmatrix(Pdim1,chi_,spliton)

if ~chi_
    mrates = zeros(Pdim1);
    if ~spliton
        % uncomment below if combined mutation and growth terms
        mrates = eye(Pdim1);
    end
else
    mrates = diag(chi_*ones(Pdim1-1,1),1) + diag(chi_*ones(Pdim1-1,1),-1);
    iloss = sum(mrates,2);
    disp(mrates);
    if ~spliton
        % uncomment below if combined mutation and growth terms
        mrates = mrates + diag(1-iloss);
    end
end
if (sum(sum(mrates<0))||(sum(iloss>1)))
    error('Qmatrix not normalized!')
end

end