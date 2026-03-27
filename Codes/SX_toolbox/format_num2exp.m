function str = format_num2exp(n)

if ischar(n)
    str = n;
else
    exp_ = floor(log10(n));
    str = sprintf('%.0fe%.0f', n/10.^exp_, exp_);
    if n==0, str='0'; end
end
end