clc
close all

SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

sigma_neuron = 20; % in deg, the width of activity distribution given each neuron (related to internal noise?)
x = -90:90;
ntrials = 1e3*4;
muORI_all = -90:1:90; % in deg, the center of (a group of) neuron(s)
nORI = length(muORI_all);
p0 = nan(ntrials, nORI);

figure('Position', [0 800 1200 200])

for iORI = 1:nORI
    muORI = muORI_all(iORI);
    sampled_activity = randn(ntrials, 1) * sigma_neuron + muORI;
    p0(:, iORI) = SX_normPDF(sampled_activity, 0, 1); % the likelihood of selecting ORI=0 given a normal dist centered at the sampled activity
    
%     subplot(1, nmu, imu), hold on
%     plot(x, SX_normPDF(x, mu_neuron, sigma_neuron), 'b-')
%     xline(0, 'k-');
%     xline(mu_neuron, 'b-');
%     if ntrials == 1, xline(sampled_activity(it), 'r-'); end
%     xlim([-90, 90])
end
% log(p0)
p0_ave = mean(p0)/sum(mean(p0));
plot(muORI_all, p0_ave)
set(findall(gcf, '-property', 'FontSize'), 'FontSize',18)
set(findall(gcf, '-property', 'LineWidth'), 'LineWidth',2)

