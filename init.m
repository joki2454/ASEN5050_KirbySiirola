%% Authors: Joshua Kirby and Amanda Siirola
%  Created: 10/26/2018
% Modified: 11/15/2018
%
% Purpose: Execute the code for the ASEN 5050 Project.
%
% Nomenclature:
%   TOF   - time of flight
%   SOI   - sphere of influence
%   CASS  - Cassini
%   CONST - constant
%   SMA   - semi-major axis
%   ECC   - eccentricity
%   INC   - inclination
%   AOP   - argument of periapsis
%   RAAN  - right ascension of the ascending node
%   nu    - true anomaly
%   T0    - initial time
%   JSOI1 - time when spacecraft enter's Jupiter's SOI
%   PJ    - time when spacecraft is at perijove
%   JSOI2 - time when spacecraft exit's Jupiter's SOI
%   SSOI  - time when spacecraft enter's Saturn's SOI
%   PS    - time when spacecraft is at perisaturnium
%   KE    - keplerian elements
%
%% Functionality Parameters
runAll     = 0;     % - run code for the entire project (necessary at least once 
                    %   if no .mat files have been previously generated and saved)
                    % - this overrides all other functionality parameters
                    
runRanging = 0;     % run code for determining allowable range for inclination

runVary    = 0;     % run code for determining DeltaV at perijove to obtain SMA 
                    % and INC in desired and determined allowable ranges,
                    % respectively
                    
runPresent = 1;     % run code for presenting results graphically

%% Housekeeping
% cleanup
close all;clc
% add paths
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
[nominal.KE_PARK,nominal.TOF,nominal.T0,nominal.JSOI1,nominal.PJ,...
  nominal.JSOI2,nominal.SSOI,nominal.PS,~] = transferSequence(nominal.DVpJ,SET);

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
presentResults(nominal,DVpJ,SET);

%% Clear all loaded SPICE kernels
% NOTE: Comment this line if you wish to use loaded MICE kernels in the command line
cspice_kclear;





















