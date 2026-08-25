RF_Track;
mass = 0.51099895069;
Q = -1;
pop = 25 * RF_Track.pC;
sigma_z = zeros(1,1000);   % preallocate
emitt_z = zeros(1,1000);   % preallocate
cut = 3;
for i = 1:1000
    % M = load(sprintf('Results/ps/phase_space_%d.dat', i));
    raw = load(sprintf('Results/ps/phase_space_%d.dat', i));
    M = reshape(raw, 6, []).';
    P = M(:,6);
    E = sqrt(P.^2 + mass^2) * 1e6;
    mE = mean(E);
    gamma = mE / mass;
    beta = sqrt(1. - 1./gamma / gamma);
    t = M(:,5) * 1e-03 * beta; % in m almost
    mP = mean(P);
    delta = (P - mP) / mP;

    mt = mean(t);
    st = std(t);
    md = mean(delta);
    sd = std(delta);
    mask = (abs(t - mt) < cut*st) & (abs(delta - md) < cut*sd);
    Z     = [t(mask) delta(mask)];
    detZ     = det(cov(Z));
    emitt_z(i)        = sqrt(abs(detZ));
end
x = 1:1000;
XSuite_gemittz_a = load('Results/XSuite_analytical_25pC_gemitt_z_1000.dat');
XSuite_gemittz_m = load('Results/XSuite_25pC_gemitt_z_1000.dat');

figure(1)
plot(x, emitt_z,'r' , 'linewidth', 2.0)
hold on
plot(x, XSuite_gemittz_m, 'b', 'linewidth', 2.0)
hold on
plot(x, XSuite_gemittz_a, 'k', 'linewidth', 2.0)
legend('RF-Track', 'Zampetakis et al.', 'Nagaitsev', 'FontSize', 20, 'location', 'northwest')
grid on
set(gca, 'FontSize', 20)
xlabel('Number of turns', 'FontSize', 20)
ylabel('Longitudinal emittance [m]', 'FontSize', 20)