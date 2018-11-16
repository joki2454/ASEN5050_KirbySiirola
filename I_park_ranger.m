%%  Author: Joshua Kirby
%  Created: 11/14/2018
% Modified: 11/15/2018
%
% Purpose: This function attempts to find solutions for DeltaV at perijove
% which produce the absolute minimum inclination (INC = 0 deg) and absolute
% maximum inclination (INC = 180 deg) using the AI_park_targeter.  This is
% done for the minimum desired SMA, maximum desired SMA, and mean desired
% SMA.  Fsolve attempts to find a solution but will always time-out
% because, with only an antivelocity impulsive maneuver at perisaturnium to
% enter parking orbits, transferSequence does not have full control over
% parking orbit inclination.  Therefore the inclination residual from the
% timed-out solutions from fsolve are considered to give bounds on the
% obtainable inclinations over the range of desired SMAs.  The ceiling rounded
% maximum of the residuals for INC = 0 deg is considered a reasonable lower
% bound on obtainable inclination.  The 180 deg plus the floor rounded minimum 
% of the residuals for  INC = 180 deg is considered a reasonable upper bound 
% on obtainnable inclination.
%
% Inputs:
%   SET       - struct of settings, initial conditions, and options
%   
% Outputs:
%   inc_range - obtainable range of inclinations when using only a single 
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
mininc = ceil(max(incres));

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
maxinc = floor(min(180 + incres));

%% Assemble range for INC
inc_range = [mininc maxinc];

end