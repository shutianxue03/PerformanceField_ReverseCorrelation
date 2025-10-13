function [x_itp,y_itp] = fxn_itp(x, y, n, flag_plot)
x_itp = linspace(min(x), max(x), n);
y_itp = interp1(x, y, x_itp, 'spline');

if flag_plot
    figure, hold on
    plot(x,y, 'o')
    plot(x_itp, y_itp, '-')
end
end