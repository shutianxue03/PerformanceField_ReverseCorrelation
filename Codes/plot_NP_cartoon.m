%% Simple sketches for neural variability and neural correlation
% Shutian-style minimal cartoons :)

clear; close all; clc

% Time axis
nT = 500;
t = linspace(-pi, 4*pi, nT);   % one cycle

% ---------- 1. Neural variability (different amplitudes) ----------
y1 = sin(t);          % neuron 1
y2 = y1/2;    % neuron 2, lower amplitude

figure('Position',[100 100 350 350]), hold on
plot(t, y1, 'k-', 'LineWidth', 2); 
plot(t, y2, 'k--', 'LineWidth', 2);
xlabel('Time')
ylabel('Activity')
ylim([-3, 3])
% title('Neuronal variability')
% axis tight
box off
set(gca, 'YTick', [], 'XTick', [])   % keep it schematic
axis off
saveas(gcf, 'NV.jpg')

%% ---------- 2. Neural correlation (identical waveform, vertical shift) ----------
offset = 1;              % vertical shift
z1 = sin(t);             % neuron 1
z2 = sin(t) - offset;    % neuron 2, same fluctuations + offset

figure('Position',[100 100 350 350]), hold on
plot(t, z1, 'k-', 'LineWidth', 2); 
plot(t, z2, 'k--', 'LineWidth', 2);
xlabel('Time')
ylabel('Activity')
% title('Neural correlation')
% axis tight
ylim([-4, 4])
box off
set(gca, 'YTick', [], 'XTick', [])
axis off

saveas(gcf, 'NC.jpg')