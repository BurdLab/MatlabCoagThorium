%
% PAPERFIGURES
%
% Matlab script to plot figures for the size-class-comparison paper
%
%
% 22-10-12: First cut:
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Adrian Burd, Dpt. of Marine Sciences, University of Georia, Athens, GA
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
clear all

x = data_select;

load rect_data.mat
load curv_data.mat

%% Plot Figures
%
% Plot figures in individual functions

% Figure 1: Particle size spectra

x = plotfigs1(diam_rect, nspec_rect, diam_curv, nspec_curv);

% Figure 2: Beta for Model A calculated from the spectral model

x = plotfigs2(t_curv, t_rect, small_2_inter_rect, small_2_inter_curv, ...
              small_2_large_rect, small_2_large_curv, mass_small_rect, mass_small_curv);

% Figure 3: betas for the Model B,

x = plotfigs3(t_rect, t_curv, beta_2_ss_curv, beta_2_ss_rect, ...
                       beta_2_sl_rect, beta_2_sl_curv);
                   
% Figure 4: Size class model solutions

x = plotfigs4(t_rect, t_curv, q_curv, q_rect, tsca_curv, tsca_rect, ...
              msca_curv, msca_rect, tscb_curv, tscb_rect, mscb_curv,...
              mscb_rect, tscc_curv, tscc_rect, mscc_curv, mscc_rect, p_rect, ...
              p_curv, vsc_curv, vsc_rect, beta_2_ss_curv, beta_2_ss_rect, ...
              beta_2_sl_rect, beta_2_sl_curv);

% Ratio of size class to sectional models          
x = plotfigs5(t_rect, t_curv, q_curv, q_rect, tsca_curv, tsca_rect, ...
              msca_curv, msca_rect, tscb_curv, tscb_rect, mscb_curv,...
              mscb_rect, tscc_curv, tscc_rect, mscc_curv, mscc_rect)

% PLot waterfall plots

x = plotfigs6(norm_data_rect, norm_data_curv);

% An additional plot

% First the rectilinear case
t_coag_b = p_rect.growth - vsc_rect(1)/p_rect.dz - beta_2_ss_rect(end)*mscb_rect(:,1);
t_coag_c = p_rect.growth - vsc_rect(1)/p_rect.dz - beta_2_ss_rect(end)*mscc_rect(:,1) - beta_2_sl_rect(end)*mscc_rect(:,2);

tg = 1/p_rect.growth;
ts = 1/(vsc_rect(1)/p_rect.dz);
tg = tg*ones(size(tscb_rect));
ts = ts*ones(size(tscb_rect));

tcb = 1./(beta_2_ss_rect(end)*mscb_rect(:,1));
tcc = 1./(beta_2_ss_rect(end)*mscc_rect(:,1) + beta_2_sl_rect(end)*mscc_rect(:,2));

figure
subplot(2,1,1)
plot(tscb_rect, tg, 'b', tscb_rect, ts, 'r', tscb_rect, tcb, 'm')
subplot(2,1,2)
plot(tscc_rect, tg, 'b', tscc_rect, ts, 'r', tscc_rect, tcc, 'm')

figure
semilogy(tscb_rect, tcb, 'm', tscc_rect, tcc, 'b--')

% Calculate slopes of spectrum

d_lionel = diam_rect(19:25);
d_lars   = diam_rect(3:26);

nr_lionel = nspec_rect(end,19:25);
nr_lars   = nspec_rect(end,3:26);

nc_lionel = nspec_curv(end,19:25);
nc_lars   = nspec_curv(end,3:26);

p_rect_lionel = polyfit(log10(d_lionel), log10(nr_lionel), 1)
p_curv_lionel = polyfit(log10(d_lionel), log10(nc_lionel), 1)
p_rect_lars   = polyfit(log10(d_lars), log10(nr_lars), 1)
p_curv_lars   = polyfit(log10(d_lars), log10(nc_lars), 1)



