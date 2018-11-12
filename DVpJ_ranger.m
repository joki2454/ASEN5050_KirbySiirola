%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   VpJ1_jc     - Velocity vector at Jupiter perijove before maneuver in
%                 Jupiter-center J2000 frame
%   SET         - Struct containing setting, initial state, and options
%   
% Outputs:
%   DV_ramrange - two-element vector containing minimum and maximum
%                 allowable DeltaV in the ram direction, km/s
%   
function [Dv_ramrange] = DVpJ_ranger(VpJ1_jc,SET)
%% Formulate parking orbit semi-major axis targeter as a function of DeltaV in the ram direction
% Designed to produce a residual between Saturn parking orbit semi-major axis and Saturn SOI
SMA_SSOItargeter_residual_eqn = @(DVpJ_ram) transferSequence(DVpJ_ram.*unitvec(VpJ1_jc),SET) - SET.CONST.SSOI/2; % km

%% Solve for DVpJ_ram to have SMA of parking orbit reach SSOI
% Using nominal case as initial guess
options = optimoptions('fsolve','maxiter',50,'TolFun',1e-4,'OptimalityTolerance',1,'TolX',1e-9,'Display',SET.TRGT.displayType);
DvpJ_ram1 = fsolve(SMA_SSOItargeter_residual_eqn,0,options); % km/s

%% Solve for DVpJ_ram to reach opposing side of SSOI
DvpJ_ram2 = fsolve(SMA_SSOItargeter_residual_eqn,-30*DvpJ_ram1,options); % km/s

%% Assemble DV_ramrange
Dv_ramrange = sort([DvpJ_ram1,DvpJ_ram2]); % km/s

end