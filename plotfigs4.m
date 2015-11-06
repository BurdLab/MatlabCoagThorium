function x = plotfigs4(t_rect, t_curv, q_curv, q_rect, tsca_curv, tsca_rect, ...
              msca_curv, msca_rect, tscb_curv, tscb_rect, mscb_curv,...
              mscb_rect, tscc_curv, tscc_rect, mscc_curv, mscc_rect, p_rect, ...
              p_curv, vsc_curv, vsc_rect, beta_2_ss_curv, beta_2_ss_rect, ...
              beta_2_sl_rect, beta_2_sl_curv)
%          
% Plot size class solutions          
%
x = 1;


psc.growth_rate = p_rect.growth;      % specific growth rate [d^{-1}]
psc.beta_ss     = beta_2_ss_rect(end);          % Coagulation rate [d^{-1}]
psc.beta_sl     = beta_2_sl_rect(end);         % Coagulation rate [d^{-1}]
psc.v           = vsc_rect;                % Settling velocity [m d^{-1}]
psc.z           = p_rect.dz;          % Layer depth [m]

modelb_ss = [(psc.growth_rate - psc.v(1)/psc.z)/psc.beta_ss (psc.growth_rate - psc.v(1)/psc.z)^2/(psc.v(2)/psc.z)/psc.beta_ss]

ps_ss = (psc.v(2)/psc.z)*(psc.growth_rate - psc.v(1)/psc.z)/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);
pl_ss = (psc.growth_rate - psc.v(1)/psc.z)^2/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);


hf = figure;

hs1 = subplot(4,1,1);
hp1 = plot(t_rect, q_rect(:,1), 'k', t_rect, (q_rect(:,2) + q_rect(:,3)), 'k--');
set(hp1(1), 'LineWidth', 2)
set(hp1(2), 'LineWidth', 2)
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,1.4e-7,'a.', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Spectral')

hs2 = subplot(4,1,2);
hp2 = plot(tsca_rect, msca_rect(:,1), 'k', tsca_rect, (msca_rect(:,2)), 'k--');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,5.6e-9,'b.', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Model A')

hs3 = subplot(4,1,3);
hp3 = plot(tscb_rect, mscb_rect(:,1), 'k', tscb_rect, (mscb_rect(:,2)), 'k--');
set(hp3(1), 'LineWidth', 2)
set(hp3(2), 'LineWidth', 2)
hold on
plot(t_rect(end), modelb_ss(1), 'ks', 'Markersize', 10, 'MarkerFaceColor', 'k')
plot(t_rect(end), modelb_ss(2), 'ks', 'Markersize', 10, 'MarkerFaceColor', 'k')
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,5.6e-7,'c.', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Model B')

hs4 = subplot(4,1,4);
hp4 = plot(tscc_rect, mscc_rect(:,1), 'k', tscc_rect, (mscc_rect(:,2)), 'k--');
set(hp4(1), 'LineWidth', 2)
set(hp4(2), 'LineWidth', 2)
hold on
plot(t_rect(end), ps_ss, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
plot(t_rect(end), pl_ss, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,0.93e-7,'d.', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Model C')

set(hs1, 'Position', [0.13 0.76726 0.775 0.186415])
set(hs2, 'Position', [0.13 0.54817 0.775 0.186415])
set(hs3, 'Position', [0.13 0.32909 0.775 0.186415])
set(hs4, 'Position', [0.13 0.110   0.775 0.186415])

orient tall
print -dpdf figure4.pdf


% get(hs1, 'Position')
% get(hs2, 'Position')
% get(hs3, 'Position')
% get(hs4, 'Position')

% Now do the same thing for the curvilinear kernels

psc.growth_rate = p_rect.growth;      % specific growth rate [d^{-1}]
psc.beta_ss     = beta_2_ss_curv(end);          % Coagulation rate [d^{-1}]
psc.beta_sl     = beta_2_sl_curv(end);         % Coagulation rate [d^{-1}]
psc.v           = vsc_curv;                % Settling velocity [m d^{-1}]
psc.z           = p_rect.dz;          % Layer depth [m]

modelb_ss = [(psc.growth_rate - psc.v(1)/psc.z)/psc.beta_ss (psc.growth_rate - psc.v(1)/psc.z)^2/(psc.v(2)/psc.z)/psc.beta_ss]

ps_ss = (psc.v(2)/psc.z)*(psc.growth_rate - psc.v(1)/psc.z)/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);
pl_ss = (psc.growth_rate - psc.v(1)/psc.z)^2/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);

hf2 = figure;

hs1 = subplot(4,1,1);
hp2 = plot(t_curv, q_curv(:,1), 'k', t_curv, (q_curv(:,2) + q_curv(:,3)), 'k--');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
set(gca, 'XTickLabel', []);
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,1.87e-6,'a.', 'FontName', 'Helvetica', 'FontSize', 18)
%set(gca, 'XTickLabel', []);

hs2 = subplot(4,1,2);
hp2 = plot(tsca_curv, msca_curv(:,1), 'k', tsca_curv, (msca_curv(:,2)), 'k--');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
set(gca, 'XTickLabel', []);
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
text(10,3.73e-7,'b.', 'FontName', 'Helvetica', 'FontSize', 18)

hs3 = subplot(4,1,3);
hp2 = plot(tscb_curv, mscb_curv(:,1), 'k', tscb_curv, (mscb_curv(:,2)), 'k--');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
hold on
plot(t_curv(end), modelb_ss(1), 'ks', 'Markersize', 10, 'MarkerFaceColor', 'k')
plot(t_curv(end), modelb_ss(2), 'ks', 'Markersize', 10, 'MarkerFaceColor', 'k')
set(gca, 'XTickLabel', []);
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', []);
text(10,5.6e-6,'c.', 'FontName', 'Helvetica', 'FontSize', 18)

hs4 = subplot(4,1,4);
hp2 = plot(tscc_curv, mscc_curv(:,1), 'k', tscc_curv, (mscc_curv(:,2)), 'k--');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
hold on
plot(t_curv(end), ps_ss, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
plot(t_curv(end), pl_ss, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'k')
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
text(10,2.8e-6,'d.', 'FontName', 'Helvetica', 'FontSize', 18)


set(hs1, 'Position', [0.13 0.76726 0.775 0.186415])
set(hs2, 'Position', [0.13 0.54817 0.775 0.186415])
set(hs3, 'Position', [0.13 0.32909 0.775 0.186415])
set(hs4, 'Position', [0.13 0.110   0.775 0.186415])

orient tall
print -dpdf figure5.pdf
