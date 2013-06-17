function mrates = Qmatrix(Pdim1,chi_,c)

mrates = zeros(Pdim1,Pdim1);
isnormz = 0;

while ~isnormz
    isnormz = 1;
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/(i-j)^c)/chi_; % mutation probability
            mrates(j,i) = mrates(i,j);
        end
        iloss = sum(mrates(i,:));
        if iloss >= 1
            isnormz = 0;
        end
        mrates(i,i) = 1-iloss;
    end
    iloss = sum(mrates(1,:));
    if iloss >= 1
        isnormz = 0;
    end
    mrates(1,1) = 1-iloss;
    isnormz
end

end