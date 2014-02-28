%%% for immune project

chi_ = 2.5e-4;
n = 400;

prefactor = (2/pi)^(1/2)*chi_;
disp([chi_,prefactor]);

xes = (1:n);
xdiffs = zeros(n);
for i=xes
    for j=i+1:n
        xdiffs(i,j) = abs(i-j)^(-2);
    end
end
mus = prefactor.*sum(xdiffs + triu(xdiffs,1)',2);
clear xdiffs;
disp(mus);

