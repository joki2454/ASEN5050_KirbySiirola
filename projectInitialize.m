%% Authors: Joshua Kirby and Amanda Siirola
%  Created: 10/26/2018
% Modified: 10/31/2018
%
% Purpose: Return execution settings for the ASEN 5050 Project.  Should be
% passed into and out of every major function in init.m.  Configurables and
% user settings should be defined here.
%
% Inputs:
% 
% Outputs: 
%   SET - struct of project settings
function SET = projectInitialize()
%% Initialization
SET = struct;

%% Filenames
SET.FILE.Iranging   = 'MAT Files/I_Park_Ranges.mat';
SET.FILE.Vary = 'MAT Files/Exploration.mat';

%% Constants
% Equatorial Radii, from App D of Vallado
SET.CONST.RS    = 60268.0; % km
SET.CONST.RJ    = 71492.0; % km

% Masses, from App D of Vallado
SET.CONST.mE    = 5.9742e24; % kg, earth mass
SET.CONST.mJ    = 318*SET.CONST.mE; % kg, jupiter mass
SET.CONST.mS    = 95.159*SET.CONST.mE; % kg, saturn mass
SET.CONST.mSun  = 332946*SET.CONST.mE; % kg, sun mass

% Gravitational Parameters, from App D of Vallado
SET.CONST.muJ   = 1.268e8; % of Jupiter, km^3/s^2
SET.CONST.muS   = 3.794e7; % of Saturn, km^3/s^2
SET.CONST.muSun = 1.32712428e11; % of Sun, km^3/s^2

%% Cassini's Nominal Flight Parameters
% Proximity Approximate Dates
SET.CASS.startDate = '20 Aug 1999 00:00:00.000';
Jdate = '30 Dec 2000';
Sdate = '01 Jul 2004';

% Sphere of Influence Radii at Proximity Approximate Date
% Jupiter SOI
temp = mice_spkezr('Jupiter Barycenter',cspice_str2et(Jdate),'J2000','NONE','Sun');
SET.CONST.JSOI = (SET.CONST.mJ/SET.CONST.mSun)^(2/5)*norm(temp.state(1:3)); % km
% Saturn SOI
temp = mice_spkezr('Saturn Barycenter',cspice_str2et(Sdate),'J2000','NONE','Sun');
SET.CONST.SSOI = (SET.CONST.mS/SET.CONST.mSun)^(2/5)*norm(temp.state(1:3)); % km

% Initial State for Cassini (two days after Earth flyby)
temp = mice_spkezr('Cassini',cspice_str2et(SET.CASS.startDate),'J2000','NONE','Sun');
SET.CASS.R0 = temp.state(1:3); % km
SET.CASS.V0 = temp.state(4:6); % km/s

% Plug SOIs and Initial State into STK, propagate using patched conics
SET.CASS.TOF_0_JSOI1    = cspice_str2et('12 Nov 2000 15:59:50.117') - cspice_str2et(SET.CASS.startDate);
SET.CASS.TOF_JSOI2_SSOI = cspice_str2et('30 Mar 2004 15:40:33.665') - cspice_str2et('19 Feb 2001 22:13:29.865');

%% Targeter Parameters (JMK kind of messy throughout project right now)
% targeter tolerances
SET.TRGT.tol = 1e-6; % km, found to be best to leave this as the default 1e-6

% Set all targeter display types
SET.TRGT.displayType = 'none'; % 'none' no fsolve display, 'final' default fsolve display

%% Solution Space Parameters
SET.RANGES.sma = [SET.CONST.RS SET.CONST.SSOI/2]; % km;
SET.RANGES.inc = [0 150]; % deg, min/max allowable, this is not the guaranteed range
                          %                         see I_park_ranger for why...

%% Presentation Parameters
% DeltaV number steps over range
SET.PRESENT.numSteps = 30;



end