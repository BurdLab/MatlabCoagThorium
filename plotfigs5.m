function x = plotfigs5(t_rect, t_curv, q_curv, q_rect, tsca_curv, tsca_rect, ...
              msca_curv, msca_rect, tscb_curv, tscb_rect, mscb_curv,...
              mscb_rect, tscc_curv, tscc_rect, mscc_curv, mscc_rect)

x = 1;
          
hf = figure

hs1 = subplot(2,2,1);
hp = plot(t_rect(80:end), mscb_rect(80:end,1)./q_rect(80:end,1), 'k', t_rect(80:end), mscb_rect(80:end,2)./(q_rect(80:end,2) + q_rect(80:end,3)), 'k--');
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
set(gca, 'Xlim', [50 300], 'YLim', [6 16])
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
legend('Small Particles', 'Large Particles', 'Location', 'SouthEast')
text(60, 15, 'a.', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')
title('Rectilinear')

hs2 = subplot(2,2,2);
hp = plot(t_curv(80:end), mscb_curv(80:end,1)./q_curv(80:end,1), 'k', t_curv(80:end), mscb_curv(80:end,2)./(q_curv(80:end,2) + q_curv(80:end,3)), 'k--');
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
set(gca, 'Xlim', [50 300], 'YLim', [1.5 3.5])
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
%ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
legend('Small Particles', 'Large Particles', 'Location', 'SouthEast')
text(60, 3.3, 'b.', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')
title('Curvilinear')

hs3 = subplot(2,2,3);
hp = plot(t_rect(80:end), mscc_rect(80:end,1)./q_rect(80:end,1), 'k', t_rect(80:end), mscc_rect(80:end,2)./(q_rect(80:end,2) + q_rect(80:end,3)), 'k--');
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'Xlim', [50 300], 'YLim', [0.5 2.0])
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
legend('Small Particles', 'Large Particles', 'Location', 'SouthEast')
text(60, 1.85, 'c.', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')
%title('Model C')

hs4 = subplot(2,2,4);
hp = plot(t_curv(80:end), mscc_curv(80:end,1)./q_curv(80:end,1), 'k', t_curv(80:end), mscc_curv(80:end,2)./(q_curv(80:end,2) + q_curv(80:end,3)), 'k--');
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'Xlim', [50 300], 'YLim', [0 2])
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
%ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
text(60, 1.8, 'd.', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')
legend('Small Particles', 'Large Particles', 'Location', 'SouthEast')
%title('Model C')

set(hs1, 'Position', [0.13    0.53884 0.36466 0.34116])
set(hs2, 'Position', [0.57034 0.53884 0.36466 0.34116])
set(hs3, 'Position', [0.13    0.11    0.36466 0.40116])
set(hs4, 'Position', [0.57034 0.11    0.36466 0.40116])

orient tall
print -dpdf figure6.pdf

get(hs1, 'Position')
get(hs2, 'Position')
get(hs3, 'Position')
get(hs4, 'Position')

