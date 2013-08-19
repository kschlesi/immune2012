function mrates = Qmatrix(Pdim1,avgsites,nu_,tstep,spliton)

%%%% abandoned so far

    mtot = poissrnd(P.*nu_.*tstep);
    newsites = floor(exprnd(avgsites,[size(mtot,1),max(mtot)])+1);
    for i=1:size(mtot,1)
        newsites(i,mtot(i)+1:end) = zeros;
    end
    %disp(newsites);
    mrates = histc(newsites',linspace(1,Pdim1,Pdim1))';
    
    
end