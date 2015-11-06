function [tsc1, msc1] = SizeClassModel1(p, t_out, q, vsc, bsc)
%
% Calculate evolution and steady state of size class model 1
%

psc.growth_rate = p.growth;      % specific growth rate [d^{-1}]
psc.coag_rate   = bsc;           % Coagulation rate [d^{-1}]
psc.v           = vsc;           % Settling velocity [m d^{-1}]
psc.z           = p.dz;          % Layer depth [m]

y_init = [q(1,1) (q(1,2) + q(1,3))];

my_options = odeset('RelTol', 1e-6, 'MaxStep', 1e-3*abs(p.t_final-p.t_init));

[tsc1, msc1] = ode23s(@(t, x) SCModel1Eqs(t, x, psc), [p.t_init p.t_final], y_init, my_options);

% Analytical Result

msc1a = SCModel1(psc, tsc1, y_init);


% Plot results

figure
hs1 = subplot(2,1,1);
hp1 = plot(tsc1, msc1(:,1), 'r', tsc1, msc1(:,2), 'b');
set(hp1(1), 'LineWidth', 2)
set(hp1(2), 'LineWidth', 2)

%xlabel('Time [d]')
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
legend('Small Particles', 'Large particles', 'Location', 'NorthWest')
hold on
hp1b = plot(tsc1, msc1a(:,1), 'r--', tsc1, msc1a(:,2), 'b--');
set(hp1b(1), 'LineWidth', 2)
set(hp1b(2), 'LineWidth', 2)
set(hs1, 'Position', [0.13 0.55 0.775 0.38])
set(gca, 'XTickLabel', [])
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)


hs2 = subplot(2,1,2);
hp2 = plot(t_out, q(:,1), 'r', t_out, (q(:,2) + q(:,3)), 'b');
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
legend('Small Particles', 'Large particles')
set(hs2, 'Position', [0.13 0.115185 0.775 0.38])
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)

% calculate size class model steady state conditions and compare with
% sectional

beta_model = p.growth - vsc(1)/p.dz;
mass_ratio = vsc(2)/p.dz/beta_model;

disp([ ])
disp(['Model A'])
disp(['Sectional: beta = ' num2str(bsc) ', mass ratio = ', num2str(q(end,1)/(q(end,2) + q(end,3)))])
disp(['Size class: beta = ' num2str(beta_model) ', mass ratio = ', num2str(mass_ratio)])
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dqdt = SCModel1Eqs(t, q, p)

dqdt = zeros(2,1);

dqdt(1) = (p.growth_rate - p.coag_rate - p.v(1)/p.z)*q(1);
dqdt(2) = p.coag_rate*q(1) - p.v(2)/p.z*q(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function msc1a = SCModel1(p, t, y0)

msc1a = zeros(length(t), 2);

a = (p.coag_rate * y0(1))/(p.growth_rate - p.coag_rate - p.v(1)/p.z + p.v(2)/p.z);
b = p.growth_rate - p.coag_rate - p.v(1)/p.z;

[p.growth_rate p.coag_rate p.v(1)]

msc1a(:,1) = y0(1)*exp(b*t);
msc1a(:,2) = (y0(2) - a)*exp(-p.v(2)/p.z*t) + a*exp(b*t);

disp(['a = ' num2str(a) ' b = ', num2str(b)])
