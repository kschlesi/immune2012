function yhat = expfun(b,x)

yhat = b(1).*exp(-1*b(2)*x);

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