function densities = unifrndpop(dim1,mean_density,lowerbound)

Nsites = mean_density*dim1/lowerbound; % number of lowerbound-sized populations to distribute
newsites = floor(unifrnd(1,dim1+1,[Nsites 1])); % randomly choose Nsites to distribute to

% distribute lowerbound-sized population to each chosen site (can choose same site twice)
densities = zeros(dim1,1);
for i=1:Nsites
    densities(newsites(i)) = densities(newsites(i))+lowerbound;
end

% the following two figures should be the same:
%
% figure
% hist(newsites,dim1)
% 
% figure
% plot((1:1:dim1),densities)

end

