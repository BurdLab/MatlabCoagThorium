function [param, opt] = SetUpCoag()
%
% SetUpCoag obtains user options for the coagulation calculations and then
% does some housecleaning
%
% USAGE:
%   
%
% HISTORY:
%   23-04-09: First cut.
%
%
% Adrian Burd, University of Georgia, 2009
%

%% Physical parameters

param.rho_fl     = 1.0275;        % Fluid density [g cm^{-3}]
param.kvisc      = 0.01;          % Kinematic viscosity [cm^2 s^{-1}]
param.g          = 980;           % Accel. due to gravity [cm s^{-2}]
param.day_to_sec = 8.64e04;       % Seconds in a day [s d^{-1}]
param.k          = 1.3e-16;       % Boltzmann's constant [erg K^{-1}]
param.r_to_rg    = 1.36;          % Interaction radius to radius of gyration

%% Section/coagulation related parameters

param.n_sections = 40;            % Nmber of sections
param.kernel     = 'KernelBrown'; % Kernel type
param.d0         = 1e-4;         % Diameter of unit particle [cm] - default 20
param.fr_dim     = 2.33;          % Particle fractal dimension
param.n1         = 10000;            % No. particles cm^{-3} in first section (40)

%% Other input parameters - set up input at a later date

param.temp    = 20 + 273;         % Temperature [K]
param.alpha   = 1.0;              % Stickiness
param.dz      = 50;               % Layer thickness [m] (65)
param.gamma   = 0.1;              % Average shear rate [s^{-1}]
param.growth  = 0.15;              % Specific growth rate in first section [d^{-1}] [0.15]
param.growth2 = 0.5;              % Higher growth rate to go to.
param.gro_sec = 2;                % Section at which growth in aggregates starts [2]
param.num_1   = 40;               % Number of particle cm^{-3} in first section


%% Paremeters for solving equations

param.t_init  = 0.0;              % Initial time for integrations [d]
param.t_final = 450.0;             % Final time for integrations [d]
param.delta_t = 1.0;              % Time interval for output [d]

%% Size Class Information
%
% First we input the desired size-class boundaries. Then we find the 
% sectional boundaries that are closest to these (see below)

param.size_class(1) = 53e-4;              % Boundary at 53 microns
param.size_class(2) = 100e-4;             % Boundary at 100 microns

param.n_size_classes = length(param.size_class);

%% Code Runtime Options

opt.tracer = 0;                   % Integrate tracer as well [0=no, 1=yes]

%% Derived parameters

param.dvisc = param.kvisc*param.rho_fl;   % Dynamic viscosity [g cm^{-1} s^{-1}]
param.del_rho = (4.5*2.48)*param.kvisc*param.rho_fl/param.g*(param.d0/2)^(-0.83);

param.conBr = 2.0/3.0*param.k*param.temp/param.dvisc;

param.a0 = param.d0/2;
param.v0 = (pi/6) * param.d0^3;

param.v_lower = param.v0*2.^( 0 : param.n_sections-1)';
param.v_upper = 2.0*param.v_lower;

param.av_vol  = 1.5*param.v_lower;
param.dcomb   = (param.v_lower*6/pi).^(1.0/3.0);
param.dwidth  = (2^(1.0/3.0)-1)*param.dcomb;

amfrac       = (4.0/3.0*pi)^(-1.0/param.fr_dim) * param.a0^(1.0-3.0/param.fr_dim);
param.amfrac = amfrac*sqrt(0.6);
param.bmfrac = 1./param.fr_dim;

param.setcon = (2.0/9.0)*param.del_rho/param.rho_fl*param.g/param.kvisc;

%% Find sectional boundaries for size class boundaries
%
% To do this we find the sectional boundary that is closest to the desired
% size-class boundary. 

% First get a vector of sectional boundaries in terms of diameter

r_i    = param.amfrac *param.v_lower.^param.bmfrac;
diam_i = 2.0*param.r_to_rg*r_i;

% Now find the section that contains each size-class boundary. Then take
% the CLOSEST sectional boundary (it may be the upper or lower boundary for 
% that section). 

for indx = 1 : param.n_size_classes

    tmp = find(diam_i < param.size_class(indx), 1, 'last');
    tmp2 = (param.size_class(indx) - diam_i(tmp))/(diam_i(tmp+1) - diam_i(tmp));
    if (tmp2 < 0.5)
        param.section(indx) = tmp;
    else
        param.section(indx) = tmp + 1;
    end
    
end




