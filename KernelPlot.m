function x = KernelPlot
%
% KernelPlot plots out the coagulation kernels for the collisions of
% particles of a given size (ra) with all the others (rb).

x = 1;

% Define a range of particle sizes

ra = [1 10 50 100 250 500 1000 1e4]*1e-4;   % Colliding particle radii [cm]
rb = logspace(-5, 0, 500);              % Target particle radii [cm]

n_ra = length(ra);                   % Number of colliding particles
n_rb = length(rb);

% Set up the constants

p = KernelSetUp(ra(1));


% Set up storage

beta_brown      = zeros(n_rb, n_ra);
beta_shear_rect = zeros(n_rb, n_ra);
beta_ds_rect    = zeros(n_rb, n_ra);
beta_shear_curv = zeros(n_rb, n_ra);
beta_ds_curv    = zeros(n_rb, n_ra);


% Calculate the kernels

for ir = 1 : n_ra
    
    beta_brown(:, ir)      = CalcBrown(ra(ir), rb, p);
    beta_shear_rect(:, ir) = CalcShearRect(ra(ir), rb, p);
    beta_ds_rect(:, ir)    = CalcDSRect(ra(ir), rb, p);
    beta_shear_curv(:, ir) = CalcShearCurv(ra(ir), rb, p);
    beta_ds_curv(:, ir)    = CalcDSCurv(ra(ir), rb, p);

end

% Now plot the kernels

rb_mat = rb';
rb_mat = rb_mat(:, ones(n_ra, 1));

ra_str = cell(size(ra));
for i_ra = 1 : n_ra
    ra_str{i_ra} = num2str(ra(i_ra));
end

figure(1)
subplot(2,2,1)
loglog(rb_mat, beta_brown)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Brownian Kernel')

subplot(2,2,2)
loglog(rb_mat, beta_shear_rect)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Rectilinear Shear Kernel')

subplot(2,2,3)
loglog(rb_mat, beta_ds_rect)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Rectilinear Differential Sedimentation Kernel')

subplot(2,2,4)
loglog(rb', beta_brown(:,1), 'b', rb', beta_shear_rect(:,1), 'r', rb', beta_ds_rect(:,1), 'g')
legend('Brownian', 'Shear', 'Diff. Sed.')
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Kernel comparison for R_a = 50 \mum')


figure(2)
subplot(2,2,1)
loglog(rb_mat, beta_brown)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Brownian Kernel')

subplot(2,2,2)
loglog(rb_mat, beta_shear_curv)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Curvilinear Shear Kernel')

subplot(2,2,3)
loglog(rb_mat, beta_ds_curv)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Curvilinear Differential Sedimentation Kernel')

subplot(2,2,4)
loglog(rb', beta_brown(:,3), 'b', rb', beta_shear_curv(:,3), 'r', rb', beta_ds_curv(:,3), 'g')
legend('Brownian', 'Shear', 'Diff. Sed.')
xlabel('R_b [cm]')
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]')
title('Kernel comparison for R_a = 50 \mum')

figure(3)
subplot(2,1,1)
semilogx(rb_mat, beta_shear_rect./beta_shear_curv)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta_{rect}/\beta_{curv}')
title('Shear Kernel')

subplot(2,1,2)
loglog(rb_mat, beta_ds_rect./beta_ds_curv)
legend(ra_str)
xlabel('R_b [cm]')
ylabel('\beta_{rect}/\beta_{curv}')
title('Differential Sedimentation Kernel')


%% Plot for paper

hf4 = figure(4);
hs1 = subplot(2,2,1);
hp1 = loglog(rb', beta_brown(:,1), 'k', rb', beta_shear_rect(:,1), 'k--', rb', beta_ds_rect(:,1), 'k-.');
legend('Brownian', 'Shear', 'Diff. Sed.', 'Location', 'NorthWest')
%xlabel('R_b [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'XTickLabel', [])
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(hp1(1), 'LineWidth', 2.0)
set(hp1(2), 'LineWidth', 2.0)
set(hp1(3), 'LineWidth', 2.0)
set(hs1, 'XTick', [1e-5 1e-4 1e-3 1e-2 1e-1 []])
set(hs1, 'YLim', [1e-14 1e6])
set(hs1, 'YTick', [1e-14 1e-10 1e-6 1e-2 1e2 1e6]);
text(4e-1, 1e-13, 'a)', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')

%title('Kernel comparison for R_a = 1 \mum')

hs2 = subplot(2,2,2);
hp2 = loglog(rb', beta_brown(:,3), 'k', rb', beta_shear_rect(:,3), 'k--', rb', beta_ds_rect(:,3), 'k-.');
legend('Brownian', 'Shear', 'Diff. Sed.', 'Location', 'NorthWest')
%xlabel('R_b [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
%ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]', 'FontName', 'Helvetica', 'FontSize', 14)
set(gca, 'XTickLabel', [])
set(gca, 'YTickLabel', [])
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(hp2(1), 'LineWidth', 2.0)
set(hp2(2), 'LineWidth', 2.0)
set(hp2(3), 'LineWidth', 2.0)
set(hs2, 'XTick', [1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
set(hs2, 'YLim', [1e-14 1e6])
set(hs2, 'YTick', [1e-14 1e-10 1e-6 1e-2 1e2 1e6]);
set(hs2, 'YTickLabel', [])
text(4e-1, 1e-13, 'b)', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')

%title('Kernel comparison for R_a = 50 \mum')

hs3 = subplot(2,2,3);
hp1 = loglog(rb', beta_brown(:,1), 'k', rb', beta_shear_curv(:,1), 'k--', rb', beta_ds_curv(:,1), 'k-.');
legend('Brownian', 'Shear', 'Diff. Sed.', 'Location', 'NorthWest')
xlabel('R_b [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]', 'FontName', 'Helvetica', 'FontSize', 18)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(hp1(1), 'LineWidth', 2.0)
set(hp1(2), 'LineWidth', 2.0)
set(hp1(3), 'LineWidth', 2.0)
set(hs3, 'XTick', [1e-5 1e-4 1e-3 1e-2 1e-1 []])
set(hs3, 'YLim', [1e-14 1e6])
%set(hs1, 'YTick', [1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4 1e6]);
set(hs3, 'YTick', [1e-14 1e-10 1e-6 1e-2 1e2 1e6]);
text(4e-1, 1e-13, 'c)', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')

%title('Kernel comparison for R_a = 1 \mum')

hs4 = subplot(2,2,4);
hp2 = loglog(rb', beta_brown(:,3), 'k', rb', beta_shear_curv(:,3), 'k--', rb', beta_ds_curv(:,3), 'k-.');
legend('Brownian', 'Shear', 'Diff. Sed.', 'Location', 'NorthWest')
xlabel('R_b [cm]', 'FontName', 'Helvetica', 'FontSize', 18)
%ylabel('\beta(R_a, R_b) [cm^{-3} s^{-1}]', 'FontName', 'Helvetica', 'FontSize', 14)
set(gca, 'FontName', 'Helvetica', 'FontSize', 18)
set(hp2(1), 'LineWidth', 2.0)
set(hp2(2), 'LineWidth', 2.0)
set(hp2(3), 'LineWidth', 2.0)
set(hs4, 'XTick', [1e-5 1e-4 1e-3 1e-2 1e-1 1e0])
set(hs4, 'YLim', [1e-14 1e6])
%set(hs2, 'YTick', [1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4 1e6]);
set(hs4, 'YTick', [1e-14 1e-10 1e-6 1e-2 1e2 1e6]);
set(hs4, 'YTickLabel', [])
text(4e-1, 1e-13, 'd)', 'FontName', 'Helvetica', 'FontSize', 18, 'FontWeight', 'bold')

%orient landscape
orient tall
%set(hs1, 'Position', [0.13 0.11 0.3856 0.815])
%set(hs2, 'Position', [0.52534 0.11 0.3856 0.815])

get(hs1, 'Position')
get(hs2, 'Position')
get(hs3, 'Position')
get(hs4, 'Position')

set(hs1, 'Position', [0.13    0.53884 0.41566 0.34116])
set(hs2, 'Position', [0.57034 0.53884 0.41566 0.34116])
set(hs3, 'Position', [0.13    0.11    0.41566 0.40116])
set(hs4, 'Position', [0.57034 0.11    0.41566 0.40116])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = KernelSetUp(r0)
% 
% Calculate the constants in the knernels and returns a structure (p)
% containing the relevant values.
%
%   r0 = radius of monomer particle [cm]
%

% Physical constants

rho_fl     = 1.0275;        % Fluid density [g cm^{-3}]
kvisc      = 0.01;          % Kinematic viscosity [cm^2 s^{-1}]
g          = 980;           % Accel. due to gravity [cm s^{-2}]
k          = 1.3e-16;       % Boltzmann's constant [erg K^{-1}]
temp       = 20 + 273;      % Temperature [K]

d0 = 2*r0;                  % Monomer diameter [cm]

dvisc = kvisc*rho_fl;       % Dynamic viscosity [g cm^{-1} s^{-1}]

del_rho = (4.5*2.48)*kvisc*rho_fl/g*(d0/2)^(-0.83);

p.br_const = 2.0/3.0*k*temp/dvisc;   % Constant for Brownian kernel

p.setcon = (2.0/9.0)*del_rho/rho_fl*g/kvisc;  % Settling velocity const

p.r_to_rg    = 1.36;          % Interaction radius to radius of gyration

p.gamma      = 0.1;              % Average shear rate [s^{-1}]


function b = CalcBrown(r0, r, p)
%
% Calculates the value of the Brownian kernel for collision between
% particles of radius r0 and all particles in sizes r.

b = (2 + r0./r + r/r0);

b = p.br_const * b(:);

function b = CalcShearRect(r0, r, p)
%
% Calclates the value of the Rectalinear Shear kernel for collision between
% particles of radius r0 and all particles in sizes r.

rg = (r0 + r)*p.r_to_rg;
b  = 1.3*rg.*rg.*rg;
b  = p.gamma*b(:);

function b = CalcDSRect(r0, r, p)
%
% Calclates the value of the Retcalinear Differential Sedimentation kernel 
% for collision between particles of radius r0 and all particles in 
% sizes r.

rg = (r0 + r)*p.r_to_rg;

v0 = SettlingVelocity(r0, r0, p.setcon);
v  = SettlingVelocity(r, r, p.setcon);

v0 = v0*ones(size(v));

b = pi * abs(v0 - v).*rg.*rg;
b = b(:);

function b = CalcShearCurv(r0, r, ps)
%
% Calclates the value of the Curvilinear Shear kernel 
% for collision between particles of radius r0 and all particles in 
% sizes r.

rg = (r0 + r)*ps.r_to_rg;

p_ratio = [r0./r; r/r0];

p  = min(p_ratio);

p1 = 1.0 + p;
p5 = p1.*p1.*p1.*p1.*p1;

efficiency = 1.0 - (1.0 + 5.0*p + 2.5*p.*p)./p5;

b = sqrt(8.0*pi/15.0)*efficiency.*rg.*rg.*rg;

b = ps.gamma*b(:);


function b = CalcDSCurv(r0, r, p)
%
% Calclates the value of the Curvilinear Shear kernel 
% for collision between particles of radius r0 and all particles in 
% sizes r.

r_mat   = [r0*ones(size(r)); r];
r_small = min(r_mat)*p.r_to_rg;

v0 = SettlingVelocity(r0, r0, p.setcon);
v  = SettlingVelocity(r, r, p.setcon);
v0 = v0*ones(size(v));

b = 0.5*pi*abs(v0 - v).*r_small.*r_small;

function v = SettlingVelocity(r, rcons, sett_const)
%
% Settling Velocity calculates the settling velocities of particles of
% given sizes.
%
% USAGE:
%   v = SettlingVelocity(r, rcons, sett_const)
%
%   v = particle settling velocities [cm s^{-1}]
%   r = particle radii [cm]
%   rcons = radii of particle conserved volumes [cm]
%   sett_const = (2g/9eta)*(delta_rho/rho_fluid)
%
% HISTORY:
%   23-04-09: First Cut
%
% Adrian Burd, University of Georgia, 2009

v = sett_const * rcons.*rcons.*rcons./r;


