function status = Qmatrix(t,y,flag,Pdim1,chi_,nu_,spliton)

global tlast ;

if strcmp(flag,'init')
    tlast = 0;
    status = 0;
else
    if numel(t)>0
        tstep = 0.1;
        %tstep = t(end)-tlast;
        %disp([t(end) tlast tstep]);
        
        %tlast = t(end);
    end
    P = y(1:Pdim1);
    mtot = poissrnd(P.*nu_.*tstep);
    newsites = chi_.*randn(max(mtot),size(mtot,1));
    newsites = floor(newsites+(newsites>=0)); 
    for i=1:size(mtot,1)
        newsites(mtot(i)+1:end,i) = zeros;
    end
    iminj = repmat(1:Pdim1,Pdim1,1)-repmat((1:Pdim1)',1,Pdim1);
    news = zeros(Pdim1,size(P,1));
    for i=1:size(P,1)
        news(:,i) = histc(newsites(:,i),iminj(i,:));
    end
    d = max(mtot).*eye(Pdim1);
    news = news - d(:,1:size(mtot,1));

    %disp(news);
    pgain = sum(news,2);
    ploss = sum(news,1)';
    disp(pgain);
    disp(ploss);
    
    mrates = news';
    status = mrates; 
    
end

end