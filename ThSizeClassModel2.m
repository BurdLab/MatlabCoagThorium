function [tscth1, thsc1, tscth2, thsc2] = ThSizeClassModel2(t_out, q, p, p2, vsc, bsc2, bsc2a, th_sc, thd)
%
% Size class model for Thorium
%


psc.growth_rate = p.growth;      % specific growth rate [d^{-1}]
psc.beta_ss     = [bsc2(287); 0.5*(bsc2(287)+bsc2(end)); bsc2(end)];     % Coagulation rate [d^{-1}]
psc.beta_sl     = [bsc2a(287); 0.5*(bsc2a(287)+bsc2a(end)); bsc2a(end)];    % Coagulation rate [d^{-1}]
psc.v           = [vsc(287,:); 0.5*(vsc(287,:)+vsc(end,:)); vsc(end,:)];           % Settling velocity [m d^{-1}]
psc.z           = p.dz;          % Layer depth [m]

psc.kdabs       = p2.kdabs;
psc.lambda_th   = p2.lambda_th    % Th decay date [d^{-1}]
psc.uran        = p2.uran;          % Uran conc [dpm cm^{-3}]
psc.th_bar      = p2.th_bar;              % Scaling factor

psc.kabs        = p2.kabs(1);

psc.prd         = p2.prd;

th_init  = [th_sc(1,1); (th_sc(1,2) + th_sc(1,3))];
y_init   = [q(1,1); (q(1,2) + q(1,3))];
thd_init = thd(1);

init_cond = [y_init; th_init; thd_init];

my_options = odeset('RelTol', 1e-8, 'MaxStep', 1e-5*abs(p.t_final-p.t_init));

[tscth1, thsc1] = ode23s(@(t, x) ThSCModel2Eqs(t, x, psc), [p.t_init : p.delta_t : p.t_final-1], init_cond, my_options);

[tscth2, thsc2] = ode23s(@(t, x) ThSCModel2bEqs(t, x, psc), [p.t_init : p.delta_t : p.t_final-1], init_cond, my_options);

% Plot the figure: analgous figure to the particle case

figure
hs = subplot(3,1,1);
hp = plot(tscth1, thsc1(:,3), 'k', tscth1, thsc1(:,4), 'k--', tscth1, thsc1(:,5), 'k-.');
hold on
set(hp(1), 'LineWidth', 2)
set(hp(2), 'LineWidth', 2)
set(hp(3), 'LineWidth', 2, 'Markersize', 10)
set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
set(gca, 'XTickLabel', [])
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
plot([290 290], [0 2e-3], 'k--')
plot([409 409], [0 2e-3], 'k--')
ylabel('Thorium Activity [dmp cm^{-03}]', 'FontName', 'Helvetica', 'FontSize', 14)
text(20,1.75e-3,'a','FontName', 'Helvetica', 'FontSize', 14)
hl = legend('Small', 'Large', 'Dissolved', 'Location', 'NorthEast');
set(hl, 'FontName', 'Helvetica', 'FontSize', 14)
%title('Model B')

hs1 = subplot(3,1,2);
hp1 = plot(tscth2, thsc2(:,3), 'k', tscth2, thsc2(:,4), 'k--', tscth2, thsc2(:,5), 'k-.');
hold on
set(hp1(1), 'LineWidth', 2)
set(hp1(2), 'LineWidth', 2)
set(hp1(3), 'LineWidth', 2, 'Markersize', 10)
set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
set(gca, 'XTickLabel', [])
%xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Thorium Activity [dmp cm^{-03}]', 'FontName', 'Helvetica', 'FontSize', 14)
plot([290 290], [0 2e-3], 'k--')
plot([409 409], [0 2e-3], 'k--')
text(20, 1.75e-3,'b','FontName', 'Helvetica', 'FontSize', 14)
hl = legend('Small', 'Large', 'Dissolved', 'Location', 'NorthEast');
set(hl, 'FontName', 'Helvetica', 'FontSize', 14)
%title('Model C')

hs2 = subplot(3,1,3);
hp2 = plot(t_out, th_sc(:,1), 'k', t_out, (th_sc(:,2) + th_sc(:,3)), 'k--', t_out, thd, 'k-.');
hold on
set(hp2(1), 'LineWidth', 2)
set(hp2(2), 'LineWidth', 2)
set(hp2(3), 'LineWidth', 2, 'Markersize', 10)
set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14)
ylabel('Thorium Activity [dmp cm^{-03}]', 'FontName', 'Helvetica', 'FontSize', 14)
plot([290 290], [0 2e-3], 'k--')
plot([409 409], [0 2e-3], 'k--')
hl=legend('Small', 'Large', 'Dissolved', 'Location', 'NorthEast');
text(20, 1.75e-3,'c','FontName', 'Helvetica', 'FontSize', 14)
set(hl, 'FontName', 'Helvetica', 'FontSize', 14)

%title('Model C')
orient tall

set(hs2, 'Position', [0.13 0.09 0.775 0.27]);
set(hs1, 'Position', [0.13 0.389 0.775 0.27]);
set(hs, 'Position', [0.13 0.688 0.775 0.27]);

disp([' '])
disp(['Model B'])
disp(['Steady State Ratios: Small Ptcles = ', num2str(thsc1(end,3)/th_sc(end,1)), ...
       ' Large particles = ', num2str(thsc1(end,4)/(th_sc(end,2) + th_sc(end,3)))])
disp(['Qs/Ql = ' num2str(th_sc(end,1)/(th_sc(end,2) + th_sc(end,3))) ', Ps/Pl = ' num2str(thsc1(end,3)/thsc1(end,4))])
disp(['Th_d(sec)/Th_d = ' num2str(thd(end)/thsc1(end,5))])

disp([' '])
disp(['Model C'])
disp(['Steady State Ratios: Small Ptcles = ', num2str(thsc2(end,3)/th_sc(end,1)), ...
       ' Large particles = ', num2str(thsc2(end,4)/(th_sc(end,2) + th_sc(end,3)))])
disp(['Qs/Ql = ' num2str(th_sc(end,1)/(th_sc(end,2) + th_sc(end,3))) ', Ps/Pl = ' num2str(thsc2(end,3)/thsc2(end,4))])
disp(['Th_d(sec)/Th_d = ' num2str(thd(end)/thsc2(end,5))])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = ThSCModel2Eqs(t, y, p)

dydt = zeros(2,1);

if t > 300
    dydt(1) = (0.5 - p.v(3,1)/p.z)*y(1) - p.beta_ss(3)*y(1)*y(1);
    dydt(2) = p.beta_ss(3)*y(1)*y(1) - p.v(3,2)/p.z*y(2);
    
    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - p.v(3,1)/p.z*y(3) - p.beta_ss(3)*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - p.v(3,2)/p.z*y(4) + p.beta_ss(3)*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);    
    
    
elseif t < 300 && t > 290
    growth = (0.5 - 0.15)/(300 - 290)*t + (0.15*300 - 0.5*290)/(300-290);
    
    v_s = (p.v(3,1) - p.v(1,1))/(300 - 290)*t + (p.v(1,1)*300 - p.v(3,1)*290)/(300-290);
    v_l = (p.v(3,2) - p.v(1,2))/(300 - 290)*t + (p.v(1,2)*300 - p.v(3,2)*290)/(300-290);
    
    b = (p.beta_ss(3) - p.beta_ss(1))/(300 - 290)*t + (p.beta_ss(1)*300 - p.beta_ss(3)*290)/(300-290);
    
    dydt(1) = (growth - v_s/p.z)*y(1) - b*y(1)*y(1);
    dydt(2) = b*y(1)*y(1) - v_l/p.z*y(2);
    
    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - v_s/p.z*y(3) - b*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - v_l/p.z*y(4) + b*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);    
        
    
else
    dydt(1) = (0.15 - p.v(1,1)/p.z)*y(1) - p.beta_ss(1)*y(1)*y(1);
    dydt(2) = p.beta_ss(1)*y(1)*y(1) - p.v(1,2)/p.z*y(2);
    
    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - p.v(1,1)/p.z*y(3) - p.beta_ss(1)*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - p.v(1,2)/p.z*y(4) + p.beta_ss(1)*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);        
    
end   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dydt = ThSCModel2bEqs(t, y, p)

dydt = zeros(2,1);

if t > 300
    dydt(1) = (0.5 - p.v(3,1)/p.z)*y(1) - p.beta_ss(3)*y(1)*y(1) - p.beta_sl(3)*y(1)*y(2);
    dydt(2) = p.beta_ss(3)*y(1)*y(1) + p.beta_sl(3)*y(1)*y(2) - p.v(3,2)/p.z*y(2);

    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - p.v(3,1)/p.z*y(3) - p.beta_ss(3)*y(1)*y(3) - p.beta_sl(3)*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - p.v(3,2)/p.z*y(4) + p.beta_ss(3)*y(1)*y(3) + p.beta_sl(3)*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);    


elseif t < 300 && t > 290
    growth = (0.5 - 0.15)/(300 - 290)*t + (0.15*300 - 0.5*290)/(300-290);
    
    v_s = (p.v(3,1) - p.v(1,1))/(300 - 290)*t + (p.v(1,1)*300 - p.v(3,1)*290)/(300-290);
    v_l = (p.v(3,2) - p.v(1,2))/(300 - 290)*t + (p.v(1,2)*300 - p.v(3,2)*290)/(300-290);
    
    b1 = (p.beta_ss(3) - p.beta_ss(1))/(300 - 290)*t + (p.beta_ss(1)*300 - p.beta_ss(3)*290)/(300-290);
    b2 = (p.beta_sl(3) - p.beta_sl(1))/(300 - 290)*t + (p.beta_sl(1)*300 - p.beta_sl(3)*290)/(300-290);

    dydt(1) = (growth - v_s/p.z)*y(1) - b1*y(1)*y(1) - b2*y(1)*y(2);
    dydt(2) = b1*y(1)*y(1) + b2*y(1)*y(2) - v_l/p.z*y(2); 
    
    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - v_s/p.z*y(3) - b1*y(1)*y(3) - b2*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - v_l/p.z*y(4) + b1*y(1)*y(3) + b2*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);    

else
    dydt(1) = (0.15 - p.v(1,1)/p.z)*y(1) - p.beta_ss(1)*y(1)*y(1) - p.beta_sl(1)*y(1)*y(2);
    dydt(2) = p.beta_ss(1)*y(1)*y(1) + p.beta_sl(1)*y(1)*y(2) - p.v(1,2)/p.z*y(2);
    
    dydt(3) = p.kabs*y(5)*y(1) - p.kdabs*y(3) - p.lambda_th*y(3) - p.v(1,1)/p.z*y(3) - p.beta_ss(1)*y(1)*y(3) - p.beta_sl(1)*y(1)*y(3);
    dydt(4) = p.kabs*y(5)*y(2) - p.kdabs*y(4) - p.lambda_th*y(4) - p.v(1,2)/p.z*y(4) + p.beta_ss(1)*y(1)*y(3) + p.beta_sl(1)*y(1)*y(3);

    dydt(5) = p.prd - p.kabs*(y(5)*y(1) + y(5)*y(2)) + p.kdabs*(y(3) + y(4)) - p.lambda_th*y(5);        

    
end   
