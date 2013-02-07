
Pdim1 = 400;
c = 0.5;
d = 5;
N = 10000;

% mtest = randn(N,1)./(d^c);
% figure
% hist(mtest,N/5)
% 
% rtest = zeros(Pdim1,1);
% for i=1:Pdim1
%     rtest(i) = randn/((i-(Pdim1/2))^c);
% end
% figure
% plot((1:Pdim1),rtest)
% 
% 
% stest = zeros(Pdim1,N);
%     for i=1:Pdim1
%         stest(i,:) = randn(N,1)./((i-(Pdim1/2))^c);
%     end
% figure
% plot((1:Pdim1),sum(stest,2)/N)


mrates = zeros(Pdim1,Pdim1);
mrates2 = zeros(Pdim1,Pdim1);
mrates3 = zeros(Pdim1,Pdim1);
    for i=2:Pdim1
        for j=1:i-1
            mrates(i,j) = (1/Pdim1)*abs(randn/((i-j)^c));
            mrates2(i,j) = abs(randn/((i-j)^c));
            mrates3(i,j) = randn/((i-j)^c);
            mrates(j,i) = mrates(i,j);
            mrates2(j,i) = mrates2(i,j);
            mrates3(j,i) = mrates3(i,j);
        end
    end
figure    
contour(mrates)
figure
contour(mrates2)
figure 
contour(mrates3)