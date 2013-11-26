function yhat = logfun(b,x)

yhat = b(2)+b(1).*log(x);

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