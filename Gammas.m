%This version gives Gamma as a gaussian
function Affinities = Gammas(y,D,Gamma_max,b)
    [m,n]=size(D);
    SquaredDistances = zeros(size(D));
    for i = 1:m
        for j = 1:n
            SquaredDistances(i,j) = (i-y(1)).^2+(j-y(2)).^2;
        end
    end
    Affinities = Gamma_max*exp(-1*SquaredDistances/(2*b^2));
end
    