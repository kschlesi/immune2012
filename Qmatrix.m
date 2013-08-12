function mrates = Qmatrix(Pdim1,chi_)

notnormz = 1;
iminj = repmat((1:Pdim1)',1,Pdim1)-repmat(1:Pdim1,Pdim1,1)+eye(Pdim1);
while notnormz
    notnormz = 0;
    mrates = (1/(Pdim1*chi_)).*abs(randn(Pdim1))./(iminj.^2).*(triu(ones(Pdim1),1));
    mrates = mrates + triu(mrates,1)';
    iloss = 1-sum(mrates,1);   % if unsymmetrising, check this!!!!
    notnormz = sum((iloss>1));
end    
mrates = mrates + diag(iloss);

end