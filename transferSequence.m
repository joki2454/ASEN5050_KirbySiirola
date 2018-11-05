%%  Author: Joshua Kirby
%  Created: 11/03/2018
% Modified: 11/03/2018
%
% Purpose: 
%
% Inputs: 
%
% Outputs:
%
function [] = transferSequence(R0,V0,TOF_0_JSOI1,DVpJ,SET)
%% Use FG Functions to propagate to Jupiter SOI
[a0,e0,i0,Omega0,omega0,nu0] = inertial2keplerian(R0,V0,SET.muS);


















end