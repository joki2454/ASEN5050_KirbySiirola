%% Authors: Joshua Kirby and Amanda Siirola
%  Created: 10/26/2018
% Modified: 11/12/2018
%
% Purpose: Execute the code for the ASEN 5050 Project.

%
%% Functionality Parameters
runAll = 0;         % - run code for the entire project (necessary at least once 
                    %   if not .mat files have been previously generated and saved)
                    % - this overrides all other functionality parameters
runRanging = 0;     % run code for determining ranges for DeltaVs applied at perijove
runVary = 1;        % run code for determining parking orbit SMA and INC for all
                    % DeltaVs within determined ranges
runPresent = 1;     % run code for presenting results graphically

%% Housekeeping
close all;clc
addpath(genpath('Subroutines'))  % Subroutines directory
addpath(genpath('mice'))         % MICE (MATLAB SPICE Interface) directory
addpath(genpath('mice kernels')) % MICE Kernel files directory
addpath(genpath('MAT Files'))    % .mat files for saving intermediate data

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
nominal.DVpJ = [0 0 0]'; % km/s
[nominal.a_park,nominal.i_park,nominal.TOF_JSOI,nominal.TOF_JSOI2_SSOI,...
  nominal.TOF_SSOI,nominal.RpJ_jc,nominal.VpJ1_jc,nominal.optim_badness] =...
    transferSequence(nominal.DVpJ,SET);

%% Range allowable DeltaVs at Perijove (along ram and anti-ram direction only for now)
if runRanging || runAll
  [inc_range] = I_park_ranger(SET);
  save(SET.FILE.Iranging,'inc_range');
else
  if isfile(SET.FILE.Iranging)
    load(SET.FILE.Iranging);
  else
    error('Attempt was made to use saved parking inclination ranges but none exist.  Set runRanging to 1 in init.m')
  end
end

% Change desired INC range in SET to be the new calculated range
inc_range = [max([inc_range(1) SET.RANGES.inc(1)]) min([inc_range(2) SET.RANGES.inc(2)])]; % deg
SET.RANGES.inc = inc_range; % deg


%% Explore parking orbit SMA and INC by varying over the range of DeltaVs
if runVary || runAll
  DVpJ = varyDVpJ(SET);
  save(SET.FILE.Vary,'DVpJ');
else
  if isfile(SET.FILE.Vary)
    load(SET.FILE.Vary);
  else
    error('Attempt was made to use saved DVpJ''s but none exist.  Set runVary to 1 in init.m')
  end
end

%% Present Results
presentResults(DVpJ,SET);

%% Clear all loaded SPICE kernels
% NOTE: Comment this line if you wish to use loaded kernels in the command line
cspice_kclear;





















