%%  Author: Joshua Kirby
%  Created: 11/05/2018
% Modified: 11/05/2018
%
% Purpose: This function produces the 3D rotation matrix from a rotating
% (r-theta-h) to cartesian (x-y-z) frame.
% 
% Use: 
%             Rxyz = C*Rrth  and Vxyz = C*Vrth
%             Rrth = C'*Rxyz and Vrth = C'*Vxyz
%
% Inputs:
%   i     - inclination, deg
%   omega - argument of periapsis, deg
%   Omega - right ascension of the ascending node, deg
%   nu    - true anomaly, deg
%
% Outputs:
%   C     - rotation matrix, 3x3
function [C] = Crth2xyz(i,omega,Omega,nu)
%% Formulate C
t = omega+nu; % deg
C = [cosd(Omega)*cosd(t)-sind(Omega)*cosd(i)*sind(t) -cosd(Omega)*sind(t)-sind(Omega)*cosd(i)*cosd(t) sind(Omega)*sind(i);...
     sind(Omega)*cosd(t)+cosd(Omega)*cosd(i)*sind(t) -sind(Omega)*sind(t)+cosd(Omega)*cosd(i)*cosd(t) -cosd(Omega)*sind(i);...
     sind(i)*sind(t)                                 sind(i)*cosd(t)                                  cosd(i)];
   




end