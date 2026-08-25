
% Read the files
rft = load('Results/RFT_25pC_nemit_x_2500_IBS_1000.dat'); % Horizontal normalized  emittance, mm * mrad
xsuite = load('Results/XSuite_25pC_nemitt_x_1000.dat'); % Horizontal normalized emittance, m * rad
xsuite_analytical = load('Results/XSuite_analytical_25pC_nemitt_x_1000.dat'); % Horizontal normalized emittance, m * rad
rft_y = load('Results/RFT_25pC_nemit_y_2500_IBS_1000.dat'); % Horizontal normalized  emittance, mm * mrad
xsuite_y = load('Results/XSuite_25pC_nemitt_y_1000.dat'); % Horizontal normalized emittance, m * rad
xsuite_analytical_y = load('Results/XSuite_analytical_25pC_nemitt_y_1000.dat'); % Horizontal normalized emittance, m * rad

% plot
% Example plotting with fixed colors and line styles
plot(rft, 'r', 'LineWidth', 2.0);                  % RF-Track x (red, solid)
hold on;
plot(xsuite * 1e6, 'b', 'LineWidth', 2.0);        % XSuite x (red, solid)
plot(xsuite_analytical * 1e6, 'k', 'LineWidth', 2.0); % XSuite analytical x (yellow, solid)

plot(rft_y, 'b-', 'LineWidth', 2.0);             % RF-Track y (red, dotted)
plot(xsuite_y * 1e6, 'r-', 'LineWidth', 2.0);   % XSuite y (red, dotted)
plot(xsuite_analytical_y * 1e6, 'k-', 'LineWidth', 2.0); % XSuite analytical y (yellow, dotted)

% Legend with epsilons
legend({
    'RF-Track \epsilon_x', 
    'Zampetakis et al. \epsilon_x', 
    'Nagaistev \epsilon_x', 
    'RF-Track \epsilon_y', 
    'Zampetakis et al. \epsilon_y', 
    'Nagaistev \epsilon_y'
}, 'FontSize', 22, 'Location', 'northwest');

xlabel('Number of turns', 'FontSize', 22)
ylabel('Normalized emittance [mm * mrad]', 'FontSize', 22)
grid on;
set(gca, 'FontSize', 22)
title('25 pC, IBS', 'FontSize', 22)

