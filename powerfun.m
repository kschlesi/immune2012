function yhat = powerfun(b,x)

yhat = b(1).*x.^(-1*b(2));

num_right = 0;
if ~all(isfinite(yhat(:)))
    disp('something bad');
    disp(yhat);
    disp(b);
    disp(num_right);
else
    num_right = num_right+1;
end

end