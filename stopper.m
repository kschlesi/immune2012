function [value,isterminal,direction] = stopper(t,y,mu_)

value = y(1:end-1)-mu_;
isterminal = ones(size(y,1)-1,1);
direction = (-1)*ones(size(y,1)-1,1);

end