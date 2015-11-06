function x = plotfigs3(t_rect, t_curv, beta_2_ss_curv, beta_2_ss_rect, ...
                       beta_2_sl_rect, beta_2_sl_curv)
% calculate the betas (or r's) for the mass loss from the small due to
% collisions within just the small size class and those between small and
% large

x = 1;

hf = figure;
hs1 = subplot(2,1,1);

[ax, h1, h2] = plotyy(t_rect, beta_2_ss_rect, t_rect, beta_2_sl_rect);
set(h1, 'LineWidth', 2, 'Color', 'k')
set(h2, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(1), 'Ycolor', 'k')
set(ax(2), 'Ycolor', 'k')
set(get(ax(1), 'Ylabel'), 'String', '$\tilde{r_1}$', 'FontName', 'Helvetica', 'FontSize', 18, 'Interpreter', 'latex')
set(get(ax(2), 'Ylabel'), 'String', '$\tilde{r_2}$', 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(1), 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(2), 'FontName', 'Helvetica', 'FontSize', 18)
hl = legend('W1', 'W2', 'Location', 'SouthEast')
h1a = findobj(get(h1,'Children'), 'String', 'W1');
set(h1a,'String', '$\tilde{r_1}$', 'Interpreter', 'latex')
h2a = findobj(get(h1,'Children'), 'String', 'W2');
set(h2a,'String', '${r}_2$', 'Interpreter', 'latex')
%set(hl, 'Interpreter', 'latex')
text(10, 3.75e5, 'a', 'FontName', 'Helvetica', 'FontSize', 18, 'Fontweight', 'bold')

hs2 = subplot(2,1,2);

[ax, h1, h2] = plotyy(t_curv, beta_2_ss_curv, t_curv, beta_2_sl_curv);
set(h1, 'LineWidth', 2, 'Color', 'k')
set(h2, 'LineWidth', 2, 'LineStyle', '--', 'Color', 'k')
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(1), 'Ycolor', 'k')
set(ax(2), 'Ycolor', 'k')
set(get(ax(1), 'Ylabel'), 'String', '$\tilde{r_1}$', 'FontName', 'Helvetica', 'FontSize', 18, 'Interpreter', 'latex')
set(get(ax(2), 'Ylabel'), 'String', '$\tilde{r_2}$', 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(1), 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(2), 'FontName', 'Helvetica', 'FontSize', 18)
hl = legend('W1', 'W2', 'Location', 'SouthEast')
%set(hl, 'Interpreter', 'latex')
h1a = findobj(get(h1,'Children'), 'String', 'W1');
set(h1a,'String', '$\tilde{r_1}$', 'Interpreter', 'latex')
h2a = findobj(get(h1,'Children'), 'String', 'W2');
set(h2a,'String', '$\tilde{r_2}$', 'Interpreter', 'latex')
text(10, 2.344e4, 'b', 'FontName', 'Helvetica', 'FontSize', 18, 'Fontweight', 'bold')

orient tall
print -dpdf figure3.pdf

get(hs1, 'Position')
get(hs2, 'Position')