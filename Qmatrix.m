function mrates = Qmatrix(Pdim1,chi_,spliton)

iminj = repmat((1:Pdim1)',1,Pdim1)-repmat(1:Pdim1,Pdim1,1)+eye(Pdim1);

if ~chi_
    mrates = zeros(Pdim1);
    if ~spliton
        % combined mutation and growth terms
        mrates = eye(Pdim1);
    end
else
    mrates = (sqrt(2/pi)/(Pdim1*chi_))./(iminj.^2).*(1-eye(Pdim1));
    iloss = sum(mrates,1);   % if unsymmetrising, check this!!!!
    if(max(iloss)>1)
        error('Mutation matrix not normalised!');
    end
    if ~spliton
        % combined mutation and growth terms
        mrates = mrates + diag(1-iloss);
    end
end

end