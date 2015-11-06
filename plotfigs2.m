function x = plotfigs2(t_curv, t_rect, small_2_inter_rect, small_2_inter_curv, ...
                       small_2_large_rect, small_2_large_curv, mass_small_rect, mass_small_curv)

% Plot the beta for model A for both rectilinear and curvilinera kernels

x = 1;

hf = figure;
hp = plot(t_rect, (small_2_inter_rect + small_2_large_rect)./(mass_small_rect), 'k-', ...
          t_curv, (small_2_inter_curv + small_2_large_curv)./(mass_small_curv), 'k--');
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
legend('Rectilinear', 'Curvilinear')
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('r_{1} [d^{-1}]', 'FontName', 'Helvetica', 'FontSize', 18)
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)

orient landscape
print -dpdf figure2.pdf
