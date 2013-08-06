
function [cabres,t] = tsnull(t,N,distrib)

% calculates average pairwise correlation and p-value
% between random vectors of m different lengths t(i)
%
% inputs: t = mx1 column vector of random vector lengths to use
%         N = number of pairs of vectors to average over, for each length
%         distrib (optional) = string specifying distribution from which 
%                              to draw random vectors (default = 'rand')
%
% outputs: cabres = mx3 matrix of average correlation (column 1),
%                   average absolute value of correlation (column 2), and
%                   average p-value (column 3) for each of the m lengths
%          t (optional) = original vector of m lengths, sorted in
%                         ascending order ; size mx1
%
% figure: plot of the 3 columns of cabres against output (sorted) t
%
% NOTE: runtime is approx. 4.5*10^-4 x N x m seconds

if nargin<3
    distrib = 'rand';
end

t = sort(t);
disp(size(t,1));
corps = zeros(size(t,1),N,2);
    for j=1:size(t,1)
        tic
        for i=1:N
            mat1 = eval(['[' distrib '(t(j),1) ' distrib '(t(j),1)]']);
            [C,P] = corrcoef(mat1);
            corps(j,i,:) = [C(1,2),P(1,2)];
        end
        toc
    end
cabs = squeeze(mean(abs(corps),2));
cres = squeeze(mean(corps,2));
cabres = [cres(:,1) cabs];

% if size(t,1)>1
%     'RUN'
%     figure
%     plot(t,cabres);
%     xlabel('number of timesteps')
%     ylabel('mean pairwise correlation')
%     legend('no abs','abs','p-value')
% end

end