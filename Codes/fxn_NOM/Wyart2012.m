
% Wyart et al., 2021, soft threshold nonlinearity

x = linspace(0,3,1e3); % since x is energy, will always be positive
alpha_all = [20:20:100]/100;
nalpha = length(alpha_all);
legends = cell(1,nalpha);

colorProp = linspace(.1, .9, length(alpha_all));
figure, hold on
for alpha = alpha_all
    ia = find(alpha == alpha_all);
    y = x+alpha.*exp(-x/alpha); plot(x,y, '-', 'color', ones(1,3)*colorProp(ia)) % p6, right volumn
    legends{ia} = sprintf('alpha = %.2f', alpha);
end

% plot(x, 2*x, 'k-') % upper limit
plot(x, x, 'r--')   % lower limit

legends = [legends, 'y=x'];
legend(legends, 'location', 'best')

xlabel('energy')
ylabel('resp variable')

set(findall(gcf, '-property', 'FontSize'), 'FontSize',15)
set(findall(gcf, '-property', 'linewidth'), 'linewidth',2)

