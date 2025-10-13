function [nn, n10] = fxn_simplifyNumber(n)
chars = num2str(n);
n10 = length(chars)-1;
nn = n/10^n10;