%%  Author: Joshua Kirby
%  Created: 10/29/2018
% Modified: 10/30/2018
%
% Purpose: 
%
% Assumptions: This function calculates the DeltaV necessary to go from a
% transfer orbit about some central body to a parking orbit about some
% other central body.  
%
% Inputs (must be SI): 
%   ap     - parking orbit semi-major axis wrt strCBp, km
%   ip     - parking orbit inclination wrt strCBp, deg
%   KEt    - 5 element vector of keplerian elements of transfer orbit wrt
%                 strCBt
%             KEt(1) = at     = semi-major axis, km
%             KEt(2) = et     = eccentricity
%             KEt(3) = it     = inclination, deg
%             KEt(4) = omegat = argument of periapsis, deg
%             KEt(5) = Omegat = right ascension of ascending node, deg
%   strCBp - string defining central body of the parking orbit
%   strCBt - string defining central body of the transfer orbit
%
% Outputs:
%
function [] = patchedConics_arrival(ap,ip,KEt,strCBp,strCBt)
%%




end