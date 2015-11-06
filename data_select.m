function x = data_select
%
% This function just loads in hte large data files for rectilinear and
% curvilinear cases, extracts the required information and then gives them
% new, unique names, saves them and clears the data space



load datafilerect.mat

t_rect = t_out;

nspec_rect = nspec;
diam_rect  = diam;
flux_rect  = flux;
mass_rect  = mass;

F_rect = F;
mass_to_rect = mass_to;

dqdt_rect = dqdt;
q_rect    = q;
vsc_rect  = vsc;
bsc_rect  = bsc;

size_class_transfer_rect = size_class_transfer;
beta_1_rect              = beta_1;
beta_2_ss_rect           = beta_2_ss;
beta_2_sl_rect           = beta_2_sl;

tsca_rect = tsca;
msca_rect = msca;
tscb_rect = tscb;
mscb_rect = mscb;
tscc_rect = tscc;
mscc_rect = mscc;

mass_small_rect = mass_small;
small_2_inter_rect = small_2_inter;
small_2_large_rect = small_2_large;

p_rect  = p;
p2_rect = p2;

norm_data_rect = norm_data;

save('rect_data.mat', 'nspec_rect', 'diam_rect', 'flux_rect', 'mass_rect', 'F_rect', ...
      'mass_to_rect', 'dqdt_rect', 'q_rect', 'vsc_rect', 'bsc_rect', 'size_class_transfer_rect', ...
      'beta_1_rect', 'beta_2_ss_rect', 'beta_2_sl_rect', 'tsca_rect', 'msca_rect', ...
      'tscb_rect', 'mscb_rect', 'tscc_rect', 'mscc_rect', 't_rect', 'mass_small_rect', ...
      'small_2_inter_rect', 'small_2_large_rect', 'p_rect', 'p2_rect', 'norm_data_rect');

clear all

% Now do the same for the curvilinear kernel. 

load datafilecurv.mat

t_curv = t_out;

nspec_curv = nspec;
diam_curv  = diam;
flux_curv  = flux;
mass_curv  = mass;

F_curv = F;
mass_to_curv = mass_to;

dqdt_curv = dqdt;
q_curv    = q;
vsc_curv  = vsc;
bsc_curv  = bsc;

size_class_transfer_curv = size_class_transfer;
beta_1_curv              = beta_1;
beta_2_ss_curv           = beta_2_ss;
beta_2_sl_curv           = beta_2_sl;

tsca_curv = tsca;
msca_curv = msca;
tscb_curv = tscb;
mscb_curv = mscb;
tscc_curv = tscc;
mscc_curv = mscc;

mass_small_curv = mass_small;
small_2_inter_curv = small_2_inter;
small_2_large_curv = small_2_large;

p_curv  = p;
p2_curv = p2;

norm_data_curv = norm_data;

save('curv_data.mat', 'nspec_curv', 'diam_curv', 'flux_curv', 'mass_curv', 'F_curv', ...
      'mass_to_curv', 'dqdt_curv', 'q_curv', 'vsc_curv', 'bsc_curv', 'size_class_transfer_curv', ...
      'beta_1_curv', 'beta_2_ss_curv', 'beta_2_sl_curv', 'tsca_curv', 'msca_curv', ...
      'tscb_curv', 'mscb_curv', 'tscc_curv', 'mscc_curv','t_curv', 'mass_small_curv', ...
      'small_2_inter_curv', 'small_2_large_curv', 'p_curv', 'p2_curv', 'norm_data_curv');

clear all

x = 1;
