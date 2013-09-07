function status = escape(t,y,flag,Pdim1,mu_,K_)

status = 0;

if strcmp(flag,'init')
    %disp(numel(t));
else
    if numel(y)
        P = y(1:Pdim1);
        if sum(P<mu_)==0  % escape
            status = 1;
        end
        if ~sum(P>=mu_)   % clearance
            status = 1;
        end
        if sum(P>=mu_)>=0.95*K_  % escape
            status = 1;
        end
    end
end

end