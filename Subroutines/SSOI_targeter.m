%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: Determines residual of an orbit's distance from Saturn's
% sphere-of-influence based on some time of flight after an initial epoch.
%
% Inputs:
%   R0   - 3D initial position vector, km
%   V0   - 3D initial velocity vector, km/s
%   mu   - gravitational parameter of central body, km^3/s^2
%   TOF  - time of flight, s
%   et0  - MICE formulated epoch at beginning of TOF, s
%   SET  - Struct containing constants and settings
%
% Outputs:
%   r    - total distance from saturn
%
function [residual] = SSOI_targeter(R0,V0,mu,TOF,et0,SET)
%% Solve for position after TOF
[R,~] = FGtime(R0,V0,mu,TOF); % km

%% Determine epoch after TOF
et = et0 + TOF; % s

%% Determine Saturn position after TOF
satstate = mice_spkezr('Saturn',et,'J2000','NONE','Sun');
Rsat = satstate.state(1:3); % km

%% Determine s/c distance from Saturn
r = norm(Rsat-R); % km

%% Determine residual distance from Saturn SOI
residual = abs(r-SET.CONST.SSOI); % km




















end