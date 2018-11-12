%%  Author: Joshua Kirby
%  Created: 11/05/2018
% Modified: 11/05/2018
%
% Purpose: Converts cccentric anomaly to true anomaly.
%
% Inputs:
%   e  - eccentricity
%   E  - eccentric anomaly, rad
%
% Outputs:
%   nu - true anomaly, deg
function [nu] = TAtoEA(e,E)
%% Convert
nu = 2*atan(sqrt((1+e)/(1-e))*tand(E/2));


end