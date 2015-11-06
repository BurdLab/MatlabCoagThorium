function [mass_change, mass_sc, vsc, bsc, mass_small, small_2_inter, small_2_large, thorium_sc] = SizeClasses(t_out, p, p2, y, mass_to, nspec, diam, flux, mass, thorium)
%
% SizeClasses lumps the sectional moreld output into size classes
%

[n_times, n_sections] = size(y);

% First calculate the mass (volume) in each size class

mass_small = sum(y(:, 1:p.section(1)-1), 2);
mass_inter = sum(y(:, p.section(1):p.section(2)-1), 2);
mass_large = sum(y(:, p.section(2):end), 2);

mass_sc = [mass_small mass_inter mass_large];

% Calculate the thorium in each size class

th_small = sum(thorium(:, 1:p.section(1)-1), 2);
th_inter = sum(thorium(:, p.section(1):p.section(2)-1), 2);
th_large = sum(thorium(:, p.section(2):end), 2);

thorium_sc = [th_small th_inter th_large];

% Calculate the total mass (volume) transfer between size classes

small_2_inter = zeros(n_times, 1);
small_2_large = zeros(n_times, 1);
inter_2_large = zeros(n_times, 1);

for i_time = 1 : n_times
    
    tmp = reshape(mass_to(:, :, i_time), n_sections, n_sections);
    
    small_2_inter(i_time) = sum(sum(tmp(p.section(1):p.section(2)-1, 1:p.section(1)-1)));
    small_2_large(i_time) = sum(sum(tmp(p.section(2):end, 1:p.section(1)-1)));
    inter_2_large(i_time) = sum(sum(tmp(p.section(2):end, p.section(1):p.section(2)-1)));
    
end

mass_change = [small_2_inter small_2_large inter_2_large];

% Calculate the total thorium transfer between size classes

th_per_part = thorium./y;

th_small_2_inter = zeros(n_times, 1);
th_small_2_large = zeros(n_times, 1);
th_inter_2_large = zeros(n_times, 1);

for i_time = 1 : n_times
    
    tmp = reshape(mass_to(:, :, i_time), n_sections, n_sections);
    
    th_part_ratio = th_per_part(i_time,:);
    th_part_ratio = th_part_ratio(ones(n_sections,1),:);
    
    tmp = tmp.*th_part_ratio;
    
    th_small_2_inter(i_time) = sum(sum(tmp(p.section(1):p.section(2)-1, 1:p.section(1)-1)));
    th_small_2_large(i_time) = sum(sum(tmp(p.section(2):end, 1:p.section(1)-1)));
    th_inter_2_large(i_time) = sum(sum(tmp(p.section(2):end, p.section(1):p.section(2)-1)));
    
end

th_change = [th_small_2_inter th_small_2_large th_inter_2_large];

% Plot some things

% figure
% semilogy(t_out, small_2_inter, 'b')
% hold on
% semilogy(t_out, small_2_large, 'r')
% semilogy(t_out, inter_2_large, 'c')
% ylabel('Rate of Transfer Between Size Classes')
% xlabel('Time')
% legend('Small to intermediate', 'Small to Large', 'Intermediate to Large', 'Location', 'SouthEast')
% 
% figure
% plot(t_out, small_2_large./(small_2_inter+small_2_large))
% ylabel('S2L/(S2L + S2I)')
% xlabel('Time')
% title('Proportion moving from small to large compared to leaving small')
% 
% figure
% plot(t_out, small_2_inter./mass_small, 'b')
% hold on
% plot(t_out, small_2_large./mass_small, 'r')'
% xlabel('Time')
% ylabel('\beta')
% legend('\beta_{s->i}', '\beta_{s->L}')
% 
% figure
% plot(t_out, (small_2_large + inter_2_large)./(mass_small + mass_inter));
% ylabel('(S2L + I2L)/(MS + MI) = \beta_1')
% xlabel('Time')
% ylabel('\beta (Two size classes)')
% 
% figure
% hp = plot(t_out, (small_2_inter + small_2_large)./(mass_small));
% set(hp, 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('\beta [d^{-1}]', 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)

% Estimate an average settling velocity for the two size classes

v_size_class = zeros(n_times, 2);

flux_size_class = [sum(flux(:, 1 : p.section(1)-1), 2) sum(flux(:, p.section(1):end), 2)];
v_size_class    = flux_size_class./[mass_small (mass_inter + mass_large)]/1e6;

% figure
% hp = plot(t_out, v_size_class(:,1), 'b', t_out, v_size_class(:,2), 'r');
% set(hp(1), 'LineWidth', 2)
% set(hp(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Settling Velocity [m d^{-1}]', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Small particles', 'Large particles')
% disp(['Small Size Class Velocity = ' num2str(v_size_class(end,1)) ' m/d'])
% disp(['Large Size Class Velocity = ' num2str(v_size_class(end,2)) ' m/d'])

vsc = v_size_class;
%vsc = [v_size_class(end,1) v_size_class(end,2)];
bsc = (small_2_inter(end) + small_2_large(end))./(mass_small(end));


