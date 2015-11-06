function x = plotfigs1(dr, nr, dc, nc)
%
% Plot size spectra
%
x = 1;

hf1 = figure(1);
hp1 = loglog(dr, nr(1, :), 'k-.', dr, nr(end, :), 'k', dc, nc(end,:), 'k--');
set(hp1(1), 'LineWidth', 1.5)
set(hp1(2), 'LineWidth', 1.5)
set(hp1(3), 'LineWidth', 1.5)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
legend('Initial', 'Final Rectilinear', 'Final Curvilinear', 'Location', 'SouthWest')
xlabel('Particle diameter [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Number spectrum [# cm^{-4}]', 'FontName', 'Helvetica', 'FontSize', 18)
axis tight

orient landscape
print -dpdf figure1.pdf
