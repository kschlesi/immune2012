function status = escape(t,y,flag,Pdim1,mu_)

status = 0;

if strcmp(flag,'init')
    %disp(numel(t));
else
    if numel(y)
        P = y(1:Pdim1);
        if sum(P<mu_)==0
            status = 1;
        end
        if ~sum(P)
            status = 1;
        end
    end
end

end