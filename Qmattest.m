%%Qmattest 
% shows histogram of (a) abs. values of N samples from a normal dist.
% with standard deviation 'chiprimeij2' (Real) vs.
% (b) N samples from the 'right side' (discarding draws from the left side) 
% of a Gaussian dist. with the same standard deviation (Paper)


clear

N = 10000000;
chiprimeij2 = 0.00025;

Qijreal = abs(chiprimeij2.*randn(N,1));
Qijpaper = -1.*ones(N,1);
for i=1:N
    while Qijpaper(i)<0
        Qijpaper(i) = chiprimeij2.*randn;
    end
end

figure
hist([Qijreal Qijpaper])
legend('Real','Paper')
