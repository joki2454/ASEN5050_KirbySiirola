%%  Author: Joshua Kirby
%  Created: 11/14/2018
% Modified: 11/14/2018
%
% Purpose: 
%
% Inputs:
%   DVpJ        - DeltaV vector in J2000 at perijove, km/s  
%   AI_park_targ- target [SMA INC]' in [km deg]
%   SET         - struct of settings, initial conditions, and options
%   
% Outputs:
%   AI_residual - residual [SMA INC]' in [km deg]
%   
function [AI_residual] = AI_park_targeter(DVpJ,AI_park_targ,SET)
%% Determine SMA and INC with DVpJ
[a_park,i_park,~,~,~,~,~,~,~] = transferSequence(DVpJ,SET);

%% Determine residual
AI_park = [a_park i_park]';
AI_residual = AI_park - AI_park_targ;
















end