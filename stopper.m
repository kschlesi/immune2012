function [value,isterminal,direction] = stopper(t,y,mu_)

value = zeros(size(y,1),1);
for i=1:size(y,1)
    value(i) = y(i)-mu_;
end
isterminal = ones(size(y,1),1);
direction = (-1)*ones(size(y,1),1);

% global Pdim1 ;       %% uncomment here to stop ONLY FOR P
% 
% value = zeros(Pdim1,1);
% for i=1:Pdim1
%     value(i) = y(i)-mu_;
% end
% isterminal = ones(Pdim1,1);
% direction = (-1)*ones(Pdim1,1);