function [tsc1, msc1, tsc2, msc2] = SizeClassModel2(t_out, q, p, vsc, bsc2, bsc2a)

psc.growth_rate = p.growth;      % specific growth rate [d^{-1}]
psc.beta_ss     = [bsc2(287); 0.5*(bsc2(287)+bsc2(end)); bsc2(end)];     % Coagulation rate [d^{-1}]
psc.beta_sl     = [bsc2a(287); 0.5*(bsc2a(287)+bsc2a(end)); bsc2a(end)];    % Coagulation rate [d^{-1}]
psc.v           = [vsc(287,:); 0.5*(vsc(287,:)+vsc(end,:)); vsc(end,:)];           % Settling velocity [m d^{-1}]
psc.z           = p.dz;          % Layer depth [m]

y_init = [q(1,1); (q(1,2) + q(1,3))];

my_options = odeset('RelTol', 1e-6, 'MaxStep', 1e-3*abs(p.t_final-p.t_init));

[tsc1, msc1] = ode23s(@(t, x) SCModel2Eqs(t, x, psc), [p.t_init : p.delta_t : p.t_final-1], y_init, my_options);

%modelb_ss = [(psc.growth_rate - psc.v(1)./psc.z)/psc.beta_ss (psc.growth_rate - psc.v(1)./psc.z).^2./(psc.v(2)./psc.z)./psc.beta_ss];

% 
% figure
% subplot(2,1,1)
% hp = plot(tsc1, msc1(:,1), 'b', tsc1, msc1(:,2), 'r');
% hold on 
% %plot(tsc1(end), modelb_ss(1), 'bs', 'MarkerSize', 10)
% %plot(tsc1(end), modelb_ss(2), 'rs', 'Markersize', 10)
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Particle Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large Particles')m 
% title('Model B')
% hs2 = subplot(2,1,2);
% hp2 = plot(t_out, q(:,1), 'b', t_out, (q(:,2) + q(:,3)), 'r');
% set(hp2(1), 'LineWidth', 2)
% set(hp2(2), 'LineWidth', 2)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large particles')
% set(hs2, 'Position', [0.13 0.115185 0.775 0.38])
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% 
% disp([' '])
% disp(['Model B'])
% disp(['Steady State Ratios: Small Ptcles (Ps/Qs) = ', num2str(msc1(end,1)/q(end,1)), ...
%       ' Large particles (Pl/Ql) = ', num2str(msc1(end,2)/(q(end,2) + q(end,3)))])
% disp(['Qs/Ql = ' num2str(q(end,1)/(q(end,2) + q(end,3))) ', Ps/Pl = ' num2str(msc1(end,1)/msc1(end,2))])
%   
y_init2 = [q(1,1); (q(1,2) + q(1,3))]; 
[tsc2, msc2] = ode23s(@(t, x) SCModel2bEqs(t, x, psc), [p.t_init : p.delta_t : p.t_final-1], y_init2, my_options);
% 
% ps_ss = (psc.v(2)/psc.z)*(psc.growth_rate - psc.v(1)/psc.z)/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);
% pl_ss = (psc.growth_rate - psc.v(1)/psc.z)^2/(psc.beta_sl*(psc.growth_rate - psc.v(1)/psc.z) + (psc.v(2)/psc.z)*psc.beta_ss);
% 
% figure
% subplot(2,1,1)
% hp = plot(tsc2, msc2(:,1), 'b', tsc2, msc2(:,2), 'r');
% hold on
% %plot(tsc2(end), ps_ss, 'bs', 'MarkerSize', 10)
% %plot(tsc2(end), pl_ss, 'rs', 'MarkerSize', 10)
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Particle Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large Particles')
% title('Model C')
% hs2 = subplot(2,1,2);
% hp2 = plot(t_out, q(:,1), 'b', t_out, (q(:,2) + q(:,3)), 'r');
% set(hp2(1), 'LineWidth', 2)
% set(hp2(2), 'LineWidth', 2)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Total Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large particles')
% set(hs2, 'Position', [0.13 0.115185 0.775 0.38])
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% 
% disp([' '])
% disp(['Model C'])
% disp(['Steady State Ratios: Small Ptcles = ', num2str(msc2(end,1)/q(end,1)), ...
%        ' Large particles = ', num2str(msc2(end,2)/(q(end,2) + q(end,3)))])
% disp(['Qs/Ql = ' num2str(q(end,1)/(q(end,2) + q(end,3))) ', Ps/Pl = ' num2str(msc2(end,1)/msc2(end,2))])
% 
% 
% figure
% hp = plot(tsc2(80:end), msc1(80:end,1)./q(80:end,1), 'b', tsc2(80:end), msc1(80:end,2)./(q(80:end,2) + q(80:end,3)), 'r');
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large Particles')
% title('Model B')
% 
% figure
% hp = plot(tsc2(10:end), msc2(10:end,1)./q(10:end,1), 'b', tsc2(10:end), msc2(10:end,2)./(q(10:end,2) + q(10:end,3)), 'r');
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('(Size Class Model)/(Sectional Model)', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small Particles', 'Large Particles')
% title('Model C')

% figure
% subplot(3,1,1)
% hp = plot(tsc1, msc1(:,1), 'b', tsc1, msc1(:,2), 'r');
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Particle Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% %title('Model B')
% 
% subplot(3,1,2)
% hp1 = plot(tsc2, msc2(:,1), 'b', tsc2, msc2(:,2), 'r');
% set(hp1(1), 'LineWidth', 2)
% set(hp1(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Particle Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% %title('Model C')
% 
% subplot(3,1,3)
% hp2 = plot(t_out, q(:,1), 'b', t_out, (q(:,2) + q(:,3)), 'r');
% set(hp2(1), 'LineWidth', 2)
% set(hp2(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Particle Volume', 'FontName', 'Helvetica', 'FontSize', 18)
% %title('Model C')
% orient tall
% h = get(hp, 'Position')
% h
% h = get(hp1, 'Position')
% h
% h = get(hp2, 'Position')
% h

function dydt = SCModel2Eqs(t, y, p)

dydt = zeros(2,1);

if t > 300
    dydt(1) = (0.5 - p.v(3,1)/p.z)*y(1) - p.beta_ss(3)*y(1)*y(1);
    dydt(2) = p.beta_ss(3)*y(1)*y(1) - p.v(3,2)/p.z*y(2);
elseif t < 300 && t > 290
    growth = (0.5 - 0.15)/(300 - 290)*t + (0.15*300 - 0.5*290)/(300-290);
    
    v_s = (p.v(3,1) - p.v(1,1))/(300 - 290)*t + (p.v(1,1)*300 - p.v(3,1)*290)/(300-290);
    v_l = (p.v(3,2) - p.v(1,2))/(300 - 290)*t + (p.v(1,2)*300 - p.v(3,2)*290)/(300-290);
    
    b = (p.beta_ss(3) - p.beta_ss(1))/(300 - 290)*t + (p.beta_ss(1)*300 - p.beta_ss(3)*290)/(300-290);
    
    dydt(1) = (growth - v_s/p.z)*y(1) - b*y(1)*y(1);
    dydt(2) = b*y(1)*y(1) - v_l/p.z*y(2);
    
else
    dydt(1) = (0.15 - p.v(1,1)/p.z)*y(1) - p.beta_ss(1)*y(1)*y(1);
    dydt(2) = p.beta_ss(1)*y(1)*y(1) - p.v(1,2)/p.z*y(2);
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = SCModel2bEqs(t, y, p)

dydt = zeros(2,1);

if t > 300
    dydt(1) = (0.5 - p.v(3,1)/p.z)*y(1) - p.beta_ss(3)*y(1)*y(1) - p.beta_sl(3)*y(1)*y(2);
    dydt(2) = p.beta_ss(3)*y(1)*y(1) + p.beta_sl(3)*y(1)*y(2) - p.v(3,2)/p.z*y(2);
elseif t < 300 && t > 290
    growth = (0.5 - 0.15)/(300 - 290)*t + (0.15*300 - 0.5*290)/(300-290);
    
    v_s = (p.v(3,1) - p.v(1,1))/(300 - 290)*t + (p.v(1,1)*300 - p.v(3,1)*290)/(300-290);
    v_l = (p.v(3,2) - p.v(1,2))/(300 - 290)*t + (p.v(1,2)*300 - p.v(3,2)*290)/(300-290);
    
    b1 = (p.beta_ss(3) - p.beta_ss(1))/(300 - 290)*t + (p.beta_ss(1)*300 - p.beta_ss(3)*290)/(300-290);
    b2 = (p.beta_sl(3) - p.beta_sl(1))/(300 - 290)*t + (p.beta_sl(1)*300 - p.beta_sl(3)*290)/(300-290);

    dydt(1) = (growth - v_s/p.z)*y(1) - b1*y(1)*y(1) - b2*y(1)*y(2);
    dydt(2) = b1*y(1)*y(1) + b2*y(1)*y(2) - v_l/p.z*y(2); 

else
    dydt(1) = (0.15 - p.v(1,1)/p.z)*y(1) - p.beta_ss(1)*y(1)*y(1) - p.beta_sl(1)*y(1)*y(2);
    dydt(2) = p.beta_ss(1)*y(1)*y(1) + p.beta_sl(1)*y(1)*y(2) - p.v(1,2)/p.z*y(2);
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

