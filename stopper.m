function [value,isterminal,direction] = stopper(t,y,Pdim1)

global mu_ ;

value = zeros(Pdim1);
for i=1:Pdim1
    value(i) = y(i)-mu_;
end
isterminal = ones(Pdim1);
direction = -1.*ones(Pdim1);