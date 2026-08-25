% Import RF-Track
RF_Track;

% Read the files
rft_IBS = load('Results/RFT_25pC_nemit_x_2500_IBS_1000.dat'); % Horizontal normalized  emittance, mm * mrad
rft_noIBS = load('Results/noIBS_emittx.dat'); % Horizontal normalized  emittance, mm * mrad

% Example plotting with fixed colors and line styles
plot(rft_IBS, 'r', 'LineWidth', 2.0);                  % RF-Track x (red, solid)
hold on;
plot(rft_noIBS, 'b', 'LineWidth', 2.0);        % XSuite x (red, solid)
legend('IBS', 'No IBS', 'FontSize', 22)
xlabel('Number of turns', 'FontSize', 22)
ylabel('Normalized emittance [mm * mrad]', 'FontSize', 22)
grid on;
set(gca, 'FontSize', 22)
title('25 pC, RF-Track simulation', 'FontSize', 22)

