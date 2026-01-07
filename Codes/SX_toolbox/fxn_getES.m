function cohenD = fxn_getES(x1, x2)

d = x1(:) - x2(:);
cohenD = mean(d,'omitnan') / std(d,0,'omitnan');