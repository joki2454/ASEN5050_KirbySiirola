%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   VpJ1_jc
%   SET
%   
% Outputs:
%   DV_ramrange
%   
function [Dv_ramrange] = DVpJ_ranger(VpJ1_jc,SET)
%% Formulate parking orbit semi-major axis targeter as a function of DeltaV in the ram direction
% Designed to produce a residual between Saturn parking orbit semi-major axis and Saturn SOI
SMA_SSOItargeter_residual_eqn = @(DVpJ_ram) transferSequence(DVpJ_ram.*unitvec(VpJ1_jc),SET) - SET.CONST.SSOI; % km

%% Solve for DVpJ_ram to have SMA of parking orbit reach SSOI
% Using nominal case as initial guess
options = optimoptions('fsolve','maxiter',50,'TolFun',1e-4,'OptimalityTolerance',1,'TolX',1e-9,'Display',SET.TRGT.displayType);
DVpJ_ram1 = fsolve(SMA_SSOItargeter_residual_eqn,0,options); % km/s

%% Solve for DVpJ_ram to reach opposing side of SSOI
DVpJ_ram2 = fsolve(SMA_SSOItargeter_residual_eqn,-30*DVpJ_ram1,options); % km/s

%% Assemble DV_ramrange
Dv_ramrange = sort([DVpJ_ram1,DVpJ_ram2]); % km/s

end