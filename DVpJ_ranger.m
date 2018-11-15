%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   RpJ         - Position vector at Jupiter perijove in Jupiter-centered
%                 J2000 frame
%   VpJ1_jc     - Velocity vector at Jupiter perijove before maneuver in
%                 Jupiter-centered J2000 frame
%   SET         - Struct containing setting, initial state, and options
%   
% Outputs:
%   DV_ramrange  - two-element vector containing minimum and maximum
%                  allowable DeltaV in the ram direction, km/s
%   Dv_rhatrange - two-element vector containing minimum and maximum 
%                  allowable DeltaV in the rhat direction, km/s
%   
function [Dv_ramrange,Dv_rhatrange] = DVpJ_ranger(RpJ_jc,VpJ1_jc,SET)
%% Solve for Dv range in the ram (theta_hat) direction
% Formulate parking orbit semi-major axis targeter as a function of DeltaV in the ram direction
% Designed to produce a residual between Saturn parking orbit semi-major axis and Saturn SOI
SMA_SSOItargeter_ram_eqn = @(DvpJ_ram) transferSequence(DvpJ_ram.*unitvec(VpJ1_jc),SET) - SET.RANGES.sma(2); % km

% Solve for DvpJ_ram to have SMA of parking orbit reach design point distance
% Using nominal case as initial guess
options = optimoptions('fsolve','maxiter',50,'TolFun',1e-4,'OptimalityTolerance',1,'TolX',1e-9,'Display',SET.TRGT.displayType);
DvpJ_ram1 = fsolve(SMA_SSOItargeter_ram_eqn,0,options); % km/s

% Solve for DVpJ_ram to reach opposing side of design point distance
DvpJ_ram2 = fsolve(SMA_SSOItargeter_ram_eqn,-1000*DvpJ_ram1,options); % km/s

% Assemble DV_ramrange
Dv_ramrange = sort([DvpJ_ram1,DvpJ_ram2]); % km/s

%% Solve for Dv range in the r_hat direction
SMA_SSOItargeter_rhat_eqn = @(DvpJ_rhat) transferSequence(DvpJ_rhat.*unitvec(RpJ_jc),SET) - SET.RANGES.sma(2); % km

% Solve for DvpJ_rhat to have SMA of parking orbit reach design point distance
% Using nominal case as initial guess
DvpJ_rhat1 = fsolve(SMA_SSOItargeter_rhat_eqn,0,options); % km/s

% Solve for DVpJ_ram to reach opposing side of design point distance
DvpJ_rhat2 = fsolve(SMA_SSOItargeter_rhat_eqn,-50*DvpJ_rhat1,options); % km/s

% Assemble DV_ramrange
Dv_rhatrange = sort([DvpJ_rhat1,DvpJ_rhat2]); % km/s
end