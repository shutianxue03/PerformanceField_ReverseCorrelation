
% Lu & Dosher, 2008, Appendix D
%% which model is this???

signal = .3; % Magnitude of the signal at the decision stage
sigma_ext = .2; % SD of external noise
sigma_int = .3; % SD of internal noise in signal-prs trials
N_abs = sigma_int;  % SD of internal noise in signal-abs trials
% U = 0;      % Number of hidden detectors
% nTrials = 1e3; % Number of trials in each pass

%% set params
signal_all = .3:.05:.5; 
nsignal = length(signal_all);
sigma_int_all = .1:.2: 1; 
nsigma_int = length(sigma_int_all);

%% empty containers
pA_all = nan(nsigma_int, nsignal);
pC_all = pA_all;
legends = cell(1,nsigma_int);

%% simulate 
for isigma_int_PRS = 1:nsigma_int % loop thru each internal noise level
    sigma_int = sigma_int_all(isigma_int_PRS);
    for idS = 1:nsignal % loop thru each signal cst 
        signal = signal_all(idS);
        %%%%%%%%%%%%%%%%%%%%%%%%
        [pC, pA] = getpC_pA(signal, sigma_ext, sigma_int);
        %%%%%%%%%%%%%%%%%%%%%%%%
        
        % So how do I simulate fake data??? how to generate responses??
        
        pA_all(isigma_int_PRS, idS) = pA;
        pC_all(isigma_int_PRS, idS) = pC;
    end
    legends{isigma_int_PRS} = sprintf('sigma_{int} = %.2f', sigma_int);
end

%% plot pC vs pA
figure, grid on
plot(pA_all', pC_all', 'o-')
xlim([.5, 1]), xticks(.5:.1:1), xlabel('pA')
ylim([.5, 1]), yticks(.5:.1:1), ylabel('pC')
legend(legends)
set(findall(gcf,'-property','FontSize'),'FontSize',12)

%% try fitting
% simulate raw data
signal_all = .3:.01:.5;
nsignal = length(signal_all);
sigma_int_all = .4;
pA_all = nan(length(sigma_int_all),nsignal );
pC_all = pA_all;
for isigma_int_PRS = 1:length(sigma_int_all)
    sigma_int = sigma_int_all(isigma_int_PRS);
    for idS = 1:nsignal
        signal = signal_all(idS);
        [pC, pA] = getpC_pA(signal, sigma_ext, sigma_int);
        pA_all(isigma_int_PRS, idS) = pA;
        pC_all(isigma_int_PRS, idS) = pC;
    end
end

%%
% for dev = [10, 50, 100, 1e3]
dev = 100;
pA_data = pA_all+randn(1,nsignal)/dev;
pC_data = pC_all+randn(1,nsignal)/dev;

%     figure, grid on, hold on
%     plot(pA_data', pC_data', 'ok')
%     plot(pA_all', pC_all', '-r')
%     xlim([.7, 1]), xticks(.7:.1:1), xlabel('pA')
%     ylim([.5, 1]), yticks(.5:.1:1), ylabel('pC')

% estimate N_int
fxn = @(N_int) getDiff(signal_all, sigma_ext, N_int, pA_data, pC_data);
N_int_est = fmincon(fxn, .3);

% predict pC and pA
pC_est_all = nan(1,nsignal);
pA_est_all = pC_est_all;
for idS = 1:nsignal
    signal = signal_all(idS);
    [pC_est, pA_est] = getpC_pA(signal, sigma_ext, N_int_est);
    pC_est_all(idS) = pC_est;
    pA_est_all(idS) = pA_est;
end

% plot
figure, grid on, hold on
plot(pA_data', pC_data', 'ok')
plot(pA_all', pC_all', '.-r') 
plot(pA_est_all', pC_est_all, '.-g')
xlim([.7, 1]), xticks(.7:.1:1), xlabel('pA')
ylim([.5, 1]), yticks(.5:.1:1), ylabel('pC')
legend('Data', 'Prediction')
%     title(dev)

%% helper fxn - getDiff
function diff = getDiff(dS_all, N_ext, N_int, pA_data, pC_data)
ndS_all = length(dS_all);
diff = 0;
for idS = 1:ndS_all
    [pC_est, pA_est] = getpC_pA(dS_all(idS), N_ext, N_int);
    diff = diff+sum((pA_data - pA_est).^2) + sum((pC_data - pC_est).^2);
end
end


%% helper fxn - getpC_pA
function [pC, pA] = getpC_pA(signal, sigma_ext, sigma_int_PRS)
% SX_normPDF = @(x,mu,sigma) exp(-(x-mu).^2/(2*sigma^2));

% get x
x = (-1e4:1e4)/100;
sigma_int_ABS = sigma_int_PRS;

% combine internal & external noise
sigma_IE = sqrt(sigma_ext^2+(sigma_int_PRS^2+sigma_int_ABS^2)/2);
t = signal + x * sigma_IE;
% get pC
pC = sum((sigma_IE)*normpdf(t,signal,sigma_int_PRS).*normcdf(t,0,sigma_int_ABS)/100);

sigma_ext2 = sigma_ext * sqrt(2);
sigma_int2 = sqrt(sigma_int_PRS^2+sigma_int_ABS^2);
t = signal + x * sigma_ext2;
% get pA
pA = sum((sigma_ext2)*normpdf(t,signal,sigma_ext2).*(normcdf(t,0,sigma_int2).^2+(1-normcdf(t,0,sigma_int2)).^2)/100);

end

