function mrates = Qmatrix(Pdim1,chi_,c)

mrates = zeros(Pdim1,Pdim1);

mrates(i,j) = 







%         for i=2:Pdim1
%             for j=1:i-1
%                 mrates(i,j) = (1/Pdim1)*abs(randn/(i-j)^c)/chi_; % mutation probability
%                 mrates(j,i) = mrates(i,j);
%             end
%             iloss = sum(mrates(i,:));
%             mrates(i,i) = 1-iloss;
%         end
% iloss = sum(mrates(1,:));
% mrates(1,1) = 1-iloss;


end