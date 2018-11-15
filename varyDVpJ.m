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
%   DVpJ
%   
function [DVpJ] = varyDVpJ(SET)
%% Solve for DVpJs over range of sma and inc
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),SET.PRESENT.numSteps); % sma
inc = linspace(SET.RANGES.inc(1),SET.RANGES.inc(2),SET.PRESENT.numSteps); % deg
perturb = [-1 -1 -1;1 1 1;1 -1 1;-1 1 1;1 -1 -1;1 1 -1]';
for i = 1:length(sma)
  for j = 1:length(inc)
    % Formulate fsolve equation
    AI_targeter_eqn = @(DVpJ) AI_park_targeter(DVpJ,[sma(i) inc(j)]',SET);
    
    % Solve
    options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display',SET.TRGT.displayType);
    [DVpJ(i,j,:),~,exitflag] = fsolve(AI_targeter_eqn,[0 0 0]',options);
    DV_firsttry = squeeze(DVpJ(i,j,:));
    
    % If no clean exit, incite fsolve to try a DeltaV in the opposite direction
    counter = 1;
    while exitflag<=0
      if counter > 5
        error('Fsolve in varyDVpJ has tried solving the same equation too many times.  Implement a better initial guess perturber.')
      end
      options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display','final');
      [DVpJ(i,j,:),~,exitflag] = fsolve(AI_targeter_eqn,perturb(:,counter).*squeeze(DV_firsttry),options);
      counter = counter + 1;
    end
  end
end