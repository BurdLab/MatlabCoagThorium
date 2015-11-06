%
% COAG_DRIVER is the main command-line driver for calculating particle size
% spectra. 
%
% USAGE:
%   Type coag_driver at the command line and follow the instructions
%
% HISTORY:
%   23-04-09: First Cut
%
% Adrian Burd, University of Georgia, 2009
%
%

close all
clear all

%% Set up and get user options
%

[p, opt] = SetUpCoag;

%% Calculate the sectionally integrated coagulation kernels
%

disp('Calculating kernels')

b_brown = CalcBetas(p);
b_brown.b1 = b_brown.b1*p.conBr*p.day_to_sec*p.alpha;
b_brown.b2 = b_brown.b2*p.conBr*p.day_to_sec*p.alpha;
b_brown.b3 = b_brown.b3*p.conBr*p.day_to_sec*p.alpha;
b_brown.b4 = b_brown.b4*p.conBr*p.day_to_sec*p.alpha;
b_brown.b5 = b_brown.b5*p.conBr*p.day_to_sec*p.alpha;

% b_brown.b1 = zeros(p.n_sections, p.n_sections);
% b_brown.b2 = zeros(p.n_sections, p.n_sections);
% b_brown.b3 = zeros(p.n_sections, p.n_sections);
% b_brown.b4 = zeros(p.n_sections, p.n_sections);
% b_brown.b5 = zeros(p.n_sections, p.n_sections);

p.kernel='KernelRectSh';
b_shear = CalcBetas(p);
b_shear.b1 = b_shear.b1*p.gamma*p.day_to_sec*p.alpha;
b_shear.b2 = b_shear.b2*p.gamma*p.day_to_sec*p.alpha;
b_shear.b3 = b_shear.b3*p.gamma*p.day_to_sec*p.alpha;
b_shear.b4 = b_shear.b4*p.gamma*p.day_to_sec*p.alpha;
b_shear.b5 = b_shear.b5*p.gamma*p.day_to_sec*p.alpha;
b_shear.b25 = b_shear.b25*p.gamma*p.day_to_sec*p.alpha;

% b_shear.b1 = zeros(p.n_sections, p.n_sections);
% b_shear.b2 = zeros(p.n_sections, p.n_sections);
% b_shear.b3 = zeros(p.n_sections, p.n_sections);
% b_shear.b4 = zeros(p.n_sections, p.n_sections);
% b_shear.b5 = zeros(p.n_sections, p.n_sections);
% b_shear.b25 = zeros(p.n_sections, p.n_sections);

p.kernel='KernelRectDS';
b_ds    = CalcBetas(p);
b_ds.b1 = b_ds.b1*p.setcon*p.day_to_sec*p.alpha;
b_ds.b2 = b_ds.b2*p.setcon*p.day_to_sec*p.alpha;
b_ds.b3 = b_ds.b3*p.setcon*p.day_to_sec*p.alpha;
b_ds.b4 = b_ds.b4*p.setcon*p.day_to_sec*p.alpha;
b_ds.b5 = b_ds.b5*p.setcon*p.day_to_sec*p.alpha;
b_ds.b25 = b_ds.b25*p.setcon*p.day_to_sec*p.alpha;

% b_ds.b1 = zeros(p.n_sections, p.n_sections);
% b_ds.b2 = zeros(p.n_sections, p.n_sections);
% b_ds.b3 = zeros(p.n_sections, p.n_sections);
% b_ds.b4 = zeros(p.n_sections, p.n_sections);
% b_ds.b5 = zeros(p.n_sections, p.n_sections);
% b_ds.b25 = zeros(p.n_sections, p.n_sections);

% Pack up the betas and store them in a new structure that will get passed
% to the derivative and jacobian calculation routines

p2.b1 =  b_brown.b1 + b_shear.b1 + b_ds.b1;
p2.b2 =  b_brown.b2 + b_shear.b2 + b_ds.b2;
p2.b3 =  b_brown.b3 + b_shear.b3 + b_ds.b3;
p2.b4 =  b_brown.b4 + b_shear.b4 + b_ds.b4;
p2.b5 =  b_brown.b5 + b_shear.b5 + b_ds.b5;

p2.b25 = p2.b2 - p2.b3 - p2.b4 - p2.b5;

%% Calculate the linear terms in the population balance equation
%

p2.growth    = CalcGrowth(p);
p2.sink_loss = CalcSinkingLoss(p);

p2.linear   = p2.growth - p2.sink_loss;

%% Initial Size Spectrum
%  Caclculate the initial size spectrum used for both estimation of the
%  steady state or evolving solution - put both in later versions

spec_init = CalcInitialSpec(p, p2);

%% Integrate Coagulation Equations
%  Set up for integrating over time

disp('Solving ODEs')

calcomp = 1:p.n_sections;
abs_tol = 1.0e-18;        % Absolute Tolerance baseline
rel_tol = 3.0e-14;        % Relative tolerance

at = (abs_tol * 1.5 .^ (-(calcomp-1)));

t_span = p.t_init : p.delta_t : p.t_final - 1;

ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at);
%ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at, 'Jacobian', @CalcCoagJac);
%ode_options = odeset('Jacobian', @CalcCoagJac);

%[t_out, y] = ode15s(@CalcCoagDeriv, t_span, spec_init, ode_options, p2); 
[t_out, y] = ode15s(@CalcCoagDeriv, t_span, spec_init, ode_options, p2); 


%% Evolve All
% Calculate the system for all the equations: first set up the Th constants

p2.kdabs     = 0.0068;        % Th deabsorption rate [d^{-1}]
p2.lambda_th = 0.028759968;   % Th decay date [d^{-1}]
p2.uran      = 2.39e-03;      % Uran conc [dpm cm^{-3}]
p2.th_bar    = 1.00;          % Scaling factor

p2.kabs      = 2.68e5*ones(p.n_sections, 1);

p2.prd = p2.lambda_th*p2.uran;

th_init = (1/(p.n_sections+1))*ones(p.n_sections+1,1);

spec_init_all = [spec_init; th_init]; 

ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', abs_tol);

[t_all_out, y_all] = ode15s(@CalcAllDeriv, t_span, spec_init_all, ode_options, p2);

% Separate out the Thorium

t_th = t_all_out;
thp_out = y_all(:,41:end-1)*p2.uran*p2.th_bar;
thd_out = y_all(:,end)*p2.uran*p2.th_bar;

%% Output
%

[nspec, diam, flux, mass] = PlotDiagnostics(p, p2, t_out, y);

%outflag = CoagOutput(p, p2, t_out, y);

[F, mass_to, norm_data] = PlotMassTo(p, p2, t_out, y);

% h2 = figure;
% movie(h2, F, 2)

%% Determine size class models

[dqdt, q, vsc, bsc, mass_small, small_2_inter, small_2_large, th_sc] = SizeClasses(t_out, p, p2, y, mass_to, nspec, diam, flux, mass, thp_out);

[size_class_transfer, beta_1, beta_2_ss, beta_2_sl] = SizeClassCoag(t_out, p, p2, y);
%[th_size_class_transfer, th_beta_1, th_beta_2_ss, th_beta_2_sl] = SizeClassCoagTh(t_out, p, p2, y,thp_out);

%% Size Class Models

[tsca, msca] = SizeClassModel1(p, t_out, q, vsc, bsc);

[tscb, mscb, tscc, mscc] = SizeClassModel2(t_out, q, p, vsc, beta_2_ss, beta_2_sl);

[tscth1, thsc1, tscth2, thsc2] = ThSizeClassModel2(t_out, q, p, p2, vsc, beta_2_ss, beta_2_sl, th_sc, thd_out);


%% Save data for later processing
%

save('datafilecurv.mat')

