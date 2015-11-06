function [nspec_i, diam_i, fluxsect, spec] = PlotDiagnostics(p, p2, t_out, spec)


%% Diagnostics
%  Run some diagnostics on the model output

[sec_gains, sec_losses]     = SectionalMassBalance(spec, p2);
[total_gains, total_losses] = TotalMassBalance(spec, p2);

%% Book-keeping
%  Here we calculate some things for the graphical output

total_spec = sum(spec');

axis_limits = [0, p.n_sections, min(total_spec), max(total_spec)];

n_times = length(t_out);

%% Calculate additional spectra
%  Calculate the number spectrum, mass spectrum and flux with respect to
%  particle diameters.

nspec_v    = zeros(n_times, p.n_sections);
massspec_v = nspec_v;
fluxsect   = nspec_v;
fluxspec   = nspec_v;

r_i = p.amfrac *p.av_vol.^p.bmfrac;
r_v = (0.75/pi*p.av_vol).^(1.0/3.0);

set_vel = SettlingVelocity(r_i, r_v, p.setcon);
set_vel = set_vel/100*p.day_to_sec;

diam_i = 2.0*p.r_to_rg*r_i;
diam_v = 2.0*r_v;

diam_i = diam_i';
diam_v = diam_v';

diam_i_mat = diam_i(ones(n_times, 1), :);
diam_v_mat = diam_v(ones(n_times, 1), :);

for jindx = 1 : n_times
    
    yout = spec(jindx,:);
    nspec_v(jindx,:)   = yout./(1.5*p.v_lower')./p.dwidth';
    masspec_v(jindx,:) = yout./p.dwidth';
    fluxsect(jindx,:)  = yout.*set_vel'*1e6;
    fluxspec(jindx,:)  = masspec_v(jindx,:).*set_vel'*1e6;
    
end
total_flux = sum(fluxsect,2);
total_mass = sum(spec, 2);

diaratio = (p.fr_dim/3)*diam_v_mat./diam_i_mat;
nspec_i = nspec_v.*diaratio;
masspec_i = masspec_v.*diaratio;
fluxspec_i = fluxspec.*diaratio;

% %% Simple spectral plot
% 
% hf1 = figure(1);
% hp1 = loglog(diam_i, nspec_i(1, :), 'b--', diam_i, nspec_i(end, :), 'r-');
% set(hp1(1), 'LineWidth', 1.5)
% set(hp1(2), 'LineWidth', 1.5)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
% legend('Initial Spectrum', 'Final Spectrum', 'FontName', 'Helvetica', 'FontSize', 18)
% xlabel('Particle diameter [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
% ylabel('Number spectrum [# cm^{-4}]', 'FontName', 'Helvetica', 'FontSize', 18)
% axis tight
% orient landscape
% 
% fig_h2 = figure(2);
% plot(t_out, total_gains.growth./(total_losses.sett + total_losses.coag))
% xlabel('Time [d]')
% ylabel('Gains/Losses')
% title('Total System Mass Balance')
% 
% fig_h3 = figure(3);
% plot(t_out, total_losses.coag./total_losses.sett)
% set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
% xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14)
% ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14)
% 
% fig_h4 = figure(4);
% plot(t_out, fluxsect, t_out, total_flux, '*--')
% xlabel('Time [d]')
% ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]')
