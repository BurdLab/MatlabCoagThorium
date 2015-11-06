function [size_class_transfer, beta_1, beta_2_ss, beta_2_sl] = SizeClassCoag(t_out, p, p2, y)

[n_times, n_sections] = size(y);

% Create a matrix of column vectors of mass - Each element of a column has 
% the same value

size_class_transfer = zeros(2, n_times);

for i_time = 1 : n_times
    
    vol = y(i_time, :);
    vol = vol(ones(n_sections, 1), :);
    
% First calculate mass transfer between size classes that arises from 
% coagulation entirely within the smaller size class. To do this we will
% use two size classes and lump the intermediate and large classes
% together. 

% First, collisions arising from collisions just within the end section

   term1 = p2.b4(p.section(1)-1, p.section(1)-1)*vol(p.section(1)-1, p.section(1)-1)*vol(p.section(1)-1, p.section(1)-1);
   
% Now the loss from the last section in the small size class arising from 
% collisions between the last section and all smaller sections

   tmp   = p2.b3 .* vol' .* vol;
   term2 = sum(tmp(:, p.section(1)-1)); 
   
% Now take account of the mass from smaller sections because of collision
% with the end section. Here we will have to subtract off the gain to the
% end section from these collisions.

   tmp   = p2.b5.*vol'.*vol - (p2.b2.*vol'.*vol)';
   term3 = sum(tmp(p.section(1)-1,:));
   
% Now look at the loss from the smaller size class because of collisions 
% with particles from the larger size class

   tmp   = p2.b5.*vol'.*vol;
   term4 = tmp(p.section(1):end, 1 : p.section(1)-1);
   term4 = sum(sum(term4));
   
   size_class_transfer(1, i_time) = term1 + term2 + term3;
   size_class_transfer(2, i_time) = term4;
end

% Figure 15: Loss from small size class because of collisions just within
% the small size class compared to those between small and large. 
% figure
% hl = plot(t_out, size_class_transfer(1,:)./(size_class_transfer(1,:) + size_class_transfer(2,:)), 'b', ...
%      t_out, size_class_transfer(2,:)./(size_class_transfer(1,:) + size_class_transfer(2,:)), 'r');
% set(hl(1), 'LineWidth', 2)
% set(hl(2), 'LineWidth', 2)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Time', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Proportion of total transfers', 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Collisions within small class', 'Collisions between small and large', 'FontName', 'Helvetica', 'FontSize', 18)
%title('Proportion of total collisions', 'FontName', 'Helvetica', 'FontSize', 18)

% Calculate the betas for these two processes assuming several things.
% First calculate the first order kinetics beta using just small x small ->
% large.



mass_small = sum(y(:, 1:p.section(1)-1), 2);
mass_large = sum(y(:, p.section(1):end), 2);

beta_1 = size_class_transfer(1,:)./mass_small';

% Now do second order kinetic terms

beta_2_ss = size_class_transfer(1,:)./mass_small'./mass_small';
beta_2_sl = size_class_transfer(2,:)./mass_small'./mass_large';


% calculate the betas (or r's) for the mass loss from the small due to
% collisions within just the small size class and those between small and
% large
figure
[ax, h1, h2] = plotyy(t_out, beta_2_ss, t_out, beta_2_sl);
set(h1, 'LineWidth', 2)
set(h2, 'LineWidth', 2, 'LineStyle', '--')
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
set(get(ax(1), 'Ylabel'), 'String', '\beta_{ss}', 'FontName', 'Helvetica', 'FontSize', 18)
set(get(ax(2), 'Ylabel'), 'String', '\beta_{sl}', 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(1), 'FontName', 'Helvetica', 'FontSize', 18)
set(ax(2), 'FontName', 'Helvetica', 'FontSize', 18)




%% Agg vs growth
% 
% Calculate the total inflow of mass to each section because of aggregation
% and compare it with the growth within that section

agg_incr = zeros(n_times, n_sections);
gro_incr = zeros(n_times, n_sections);

for i_time = 1 : n_times
    
    vcon_r = y(i_time, :);
    vcon   = vcon_r';
    
    vcon_shift = [0 vcon_r(1:n_sections-1)];

    term1 = vcon_r * p2.b2;
    term1 = vcon_r .* term1;

    term2 = vcon_r * p2.b1;
    term2 = term2 .* vcon_shift;

    term3 = p.growth * vcon;
    
    agg_incr(i_time, :) = (term1 + term2);
    gro_incr(i_time, :) = term3';
    
end


x  = [1 : n_sections]';
x2 = x(:, ones(1,n_times));
t  = t_out';
t2 = t(ones(n_sections,1), :);
figure
contourf(t2, x2, log10(gro_incr'./agg_incr'))
colormap('pink')
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
colorbar
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('Section Number', 'FontName', 'Helvetica', 'FontSize', 18)



