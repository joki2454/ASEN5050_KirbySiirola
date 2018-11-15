%%  Author: Joshua Kirby
%  Created: 11/14/2018
% Modified: 11/14/2018
%
% Purpose: 
%
% Inputs:
%   SET
%   
% Outputs:
%   inc_range
%   
function [inc_range] = I_park_ranger(SET)
%% Determine minimum for INC
% Options for fsolve
options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display',SET.TRGT.displayType);

% SMA's to check
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),3); % km

% Attempt to solve with minimum inclination
for i = 1:length(sma)
  % SMA INC targeter with minimum inclination (inc = 0 deg)
  AI_targeter_eqn = @(DVpJ) AI_park_targeter(DVpJ,[sma(i) 0]',SET);

  % Determine residuals
  [~,res(:,i)] = fsolve(AI_targeter_eqn,[0 0 0]',options);
end

% Determine inclination residuals
incres = res(2,:); % deg

% Identify reasonable minimum inclination
mininc = ceil(max(incres) + 1);

%% Determine maximum for INC
% Attempt to solve with maximum inclination
for i = 1:length(SET.RANGES.sma)
  % SMA INC targeter with minimum inclination (inc = 180 deg)
  AI_targeter_eqn = @(DVpJ) AI_park_targeter(DVpJ,[SET.RANGES.sma(i) 180]',SET);

  % Determine residuals
  [~,res(:,i)] = fsolve(AI_targeter_eqn,[0 0 0]',options);
end

% Determine inclination residuals
incres = res(2,:); % deg

% Identify reasonable minimum inclination
maxinc = floor(min(180 + incres - 1));

%% Assemble range for INC
inc_range = [mininc maxinc];

end