%%  Author: Joshua Kirby
%  Created: 11/14/2018
% Modified: 11/15/2018
%
% Purpose: This solves for the necessary DeltaV performed at perijove to
% obtain a desired SMA and INC for the parking orbit about Saturn.  This is
% done using fsolve and the AI_park_targeter subroutine.  A method for
% perturbing the initial condition in the case of fsolve failure is
% implemented.  The perturbation method is as follows: 1) initial guess of
% DVpJ = [0 0 0]' km/s, 2) initial guess of an adjacent solution if it
% exists, 3) initial guess of the result from 1 rotated in random
% directions until, 4) after five tries give up and throw an error.
%
% Inputs:
%   SET   - struct of settings, initial conditions, and options
%   
% Outputs:
%   DVpJ  - SET.PRESENT.numSteps x SET.PRESENT.numSteps x 3 matrix of
%           DeltaV vectors in J2000 performed at perijove to attain a
%           desired sma and inc (each corresponding to the first and second
%           dimension of the matrix, respectively), km/s
%   
function [DVpJ] = varyDVpJ(SET)
%% Allocation
% Produce sma and inc vectors
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),SET.PRESENT.numSteps); % sma
inc = linspace(SET.RANGES.inc(1),SET.RANGES.inc(2),SET.PRESENT.numSteps); % deg

% Set of random directions used to perturb the initial guess for fsolve
perturb = [-1 -1 -1;1 -1 1;-1 1 1;1 -1 -1;1 1 -1]';

% preallocate DVpJ
DVpJ = zeros(length(sma),length(inc),3);

%% Solve for DVpJs over range of sma and inc
for i = 1:length(sma)
  for j = 1:length(inc)
    % Formulate fsolve equation
    AI_targeter_eqn = @(DVpJ) AI_park_targeter(DVpJ,[sma(i) inc(j)]',SET);
    
    % Start some tracking variables
    exitflag = -1;
    counter = 1;
    
    % Attempt to solve until solved (or too many attempts)
    while exitflag <= 0 % exitflag <= 0 means fsolve failed to solve
      % First try an initial guess of DVpJ = [0 0 0]' km/s and don't print results
      if (counter == 1)
        options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display','none');
        [DVpJ(i,j,:),~,exitflag] = fsolve(AI_targeter_eqn,[0 0 0]',options);
        DV_firsttry = squeeze(DVpJ(i,j,:)); % save result of first attempt
      
        % next try using an adjacent solution (if one exists) as an initial 
      %      guess, print result if desired
      elseif counter == 2 && (j > 1) && (i > 1)
        options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display',SET.TRGT.displayType);
        [DVpJ(i,j,:),~,exitflag] = fsolve(AI_targeter_eqn,squeeze(DVpJ(i-1,j-1,:)),options);
      
        % if fsolve has tried 5 times and failed then this while loop isn't working
      elseif counter > 5
        error('Fsolve in varyDVpJ has tried solving the same equation too many times.  Implement a better initial guess perturber.')
      
        % rotate the result of the first attempt to some other direction
        %        and try again
      elseif counter > 1
        options = optimoptions('fsolve','algorithm','levenberg-marquardt','TolFun',1e-6,'Display',SET.TRGT.displayType);
        [DVpJ(i,j,:),~,exitflag] = fsolve(AI_targeter_eqn,perturb(:,counter-1).*DV_firsttry,options);
      
      % If counter is none of the above values... that's bad
      else
        error('Invalid counter value.  Fix this code.')
      end % if
      counter = counter + 1;
    end % while
  end % for inc
end % for sma
      
end