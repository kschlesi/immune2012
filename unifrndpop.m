function densities = unifrndpop(dim1,mean_density,width)

densities = unifrnd(mean_density-width/2,mean_density+width/2,[dim1 1]);
% figure
% plot((1:1:dim1),densities)

end

