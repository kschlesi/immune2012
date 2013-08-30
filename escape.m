function status = escape(t,y,flag,Pdim1,mu_)

status = 0;
if strcmp(flag,'init')
    disp(numel(t));
else
    P = y(1:Pdim1);
    if sum(P<mu_)==0
        status = 1;
    end
end

end