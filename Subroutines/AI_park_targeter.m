%%  Author: Joshua Kirby
%  Created: 11/14/2018
% Modified: 11/15/2018
%
% Purpose: This function uses transferSequence to determine the Saturn
% parking orbit SMA and INC for a given DeltaV at perijove and produces the
% residual between the obtained SMA and INC and the desired SMA and INC.
% It is used as the targeter function with which fsolve determines DVpJ.
%
% Inputs:
%   DVpJ        - DeltaV vector in J2000 at perijove, km/s  
%   AI_park_targ- target parking orbit [SMA INC]' in [km deg]
%   SET         - struct of settings, initial conditions, and options
%   
% Outputs:
%   AI_residual - residual [SMA INC]' in [km deg]
%   
function [AI_residual] = AI_park_targeter(DVpJ,AI_park_targ,SET)
%% Determine SMA and INC with DVpJ
KE_PARK = transferSequence(DVpJ,SET);

%% Determine residual
AI_park = [KE_PARK.a KE_PARK.i]';
AI_residual = AI_park - AI_park_targ;
















end