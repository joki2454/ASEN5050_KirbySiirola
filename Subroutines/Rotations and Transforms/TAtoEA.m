%%  Author: Joshua Kirby
%  Created: 11/05/2018
% Modified: 11/05/2018
%
% Purpose: Converts true anomaly to eccentric anomaly.
%
% Inputs:
%   e  - eccentricity
%   nu - true anomaly, deg
%
% Outputs:
%   E  - eccentric anomaly, rad
function [E] = TAtoEA(e,nu)
%% Convert
E = 2*atan(sqrt((1-e)/(1+e))*tand(nu/2));


end