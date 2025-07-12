function cohenD = fxn_getES(x1, x2)

m1 = mean(x1);
m2 = mean(x2);
var1 = var(x1);
var2 = var(x2);
cohenD = abs(m1-m2)/sqrt((var1+var2)/2);