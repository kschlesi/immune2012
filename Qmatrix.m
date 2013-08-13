function mrates = Qmatrix(Pdim1,chi_,w_)

% inputs: w_ is a column vector giving the fraction of pathogen that mutate
%         away from each shape space site -- a.k.a. the MUTATION LANDSCAPE.

%w_ = 6.4*10^-5*ones(Pdim1,1);
%w_ = 0.1*(1:Pdim1)';

iminj = repmat((1:Pdim1)',1,Pdim1)-repmat(1:Pdim1,Pdim1,1)+eye(Pdim1);
if chi_
    mrates = (exp(iminj/chi_)).*(triu(ones(Pdim1),1));
    mrates = mrates + triu(mrates,1)';
    mrates = mrates./repmat(sum(mrates,2)./w_,1,5); % if unsymmetrising, check!!
    %disp([w_ sum(mrates,2)]);
    % uncomment below if combined mutation and growth terms
    % mrates = mrates + diag(1-w_);
    % disp(sum(mrates,2));
else
    mrates = zeros(Pdim1);
    % uncomment below if combined mutation and growth terms
    % mrates = eye(Pdim1);
end

end