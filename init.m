%% Authors: Joshua Kirby and Amanda Siirola
%  Created: 10/26/2018
% Modified: 10/30/2018
%
% Purpose: Execute the code for the ASEN 5050 Project.
%
%% Housekeeping
close all;clc
addpath(genpath('Subroutines'))  % Subroutines directory
addpath(genpath('mice'))         % MICE (MATLAB SPICE Interface) directory
addpath(genpath('micekernels')) % MICE Kernel files directory

%% Initialize MICE (MATLAB SPICE Interface) Kernels
kernelpath = './mice kernels/Cassini Kernels/'; % directory for Cassini Kernels

filestrs = {'cas_1999_v22_990816_990820.tm',...
            'cas_2000_v24_001228_010101.tm',...
            'cas_2001_v25_001228_010101.tm',...
            'cas_2004_v26_040629_040703.tm'}; % Cassini Kernels
for i = 1:length(filestrs)
  filestrs{i} = [kernelpath filestrs{i}];
end

% load kernels
cspice_furnsh(filestrs); % Cassini kernels
cspice_furnsh('./mice kernels/Planetary Ephemeris/de432s.bsp'); % planetary kernel

%% Project Settings
SET = projectInitialize();

%% Define Nominal Orbit (No maneuvers at perijove)


%% Explore SMA and I arrival solution space (ignore RAAN? that's what I'm thinking)


%% Present Results


%% Clear all loaded SPICE kernels
% NOTE: Stop before this line if you wish to use loaded kernels in the
%       command line
%cspice_kclear;