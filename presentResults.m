%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/15/2018
%
% Purpose:  This function presents various results for the project in the
% form of plots.
%
% Inputs:
%   nominal - struct containing all data from the nominal (DVpJ = [0 0 0]'
%             km/s) case
%   DVpJ    - SET.PRESENT.numSteps x SET.PRESENT.numSteps x 3 matrix of
%             DeltaV vectors in J2000 performed at perijove to attain a
%             desired sma and inc (each corresponding to the first and second
%             dimension of the matrix, respectively), km/s
%   SET     - struct of settings, initial conditions, and options
%   
function [] = presentResults(nominal,DVpJ,SET)
%% Calculate sma and inc vectors and produce meshgrid
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),SET.PRESENT.numSteps); % sma
inc = linspace(SET.RANGES.inc(1),SET.RANGES.inc(2),SET.PRESENT.numSteps); % deg
[SMA,INC] = meshgrid(sma,inc);

%% Calculate transfer parameters over entire range
for i = 1:size(DVpJ,1)
  for j = 1:size(DVpJ,2)
    [KE_PARK(i,j),TOF(i,j),T0(i,j),JSOI1(i,j),PJ(i,j),JSOI2(i,j),SSOI(i,j),PS(i,j),optim_badness] = ...
          transferSequence(squeeze(DVpJ(i,j,:)),SET);
    % Make sure all solutions actually enter Saturn's SOI, otherwise there
    %   is a serious problem
    if optim_badness
      error(['One of the solutions for DVpJ produces an orbit which does not ',...
              'intersect Saturn''s sphere-of-influence, which is not handled by ',...
              'the code as-is.  This most likely occurred because the upper range ',...
              'for target semi-major axes (SET.RANGES.sma(2) in projectInitialize.m) ',...
              'is too large, try reducing it to something ',...
              'significantly less than Saturn''s sphere-of-influence'])
    end
  end
end

%% Calculate all |DeltaV|s and residuals
for i = 1:length(sma)
  for j = 1:length(inc)
    res(i,j)    = sqrt(((KE_PARK(i,j).a-sma(i))/sma(i))^2 + ((KE_PARK(i,j).i-inc(j))/inc(j))^2); % %
    Dv_pJ(i,j)  = norm(squeeze(DVpJ(i,j,:)));                    % km/s
    Dv_pS(i,j)  = norm(PS(i,j).DV);                              % km/s
    Dv_tot(i,j) = norm(squeeze(DVpJ(i,j,:))) + norm(PS(i,j).DV); % km/s
  end
end

%% |DeltaV| at perijove
figure('name','Dv_pJ')
hold on
grid on
grid minor
surf(SMA,INC,Dv_pJ)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{perijove} (km/s)')
title('|\DeltaV|_{pJ} vs. SMA and INC')
hold off

%% |DeltaV| at perisaturnium
figure('name','Dv_pS')
hold on
grid on
grid minor
surf(SMA,INC,Dv_pS)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{perisaturnium} (km/s)')
title('|\DeltaV|_{pS} vs. SMA and INC')
hold off

%% |DeltaV| total
figure('name','Dv_tot')
hold on
grid on
grid minor
surf(SMA,INC,Dv_tot)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('|\DeltaV|_{total} (km/s)')
title('|\DeltaV|_{tot} vs. SMA and INC')
hold off

%% Sum of squared residuals
figure('name','residuals')
hold on
grid on
grid minor
surf(SMA,INC,res)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('(\deltaSMA^2+\deltaINC^2)^(1/2) []')
title('Sum of Squared SMA/INC Residuals vs. SMA and INC')
hold off

%% RAAN
figure('name','RAAN')
hold on
grid on
grid minor
for i = 1:size(KE_PARK,2)
  raan(:,i) = extractfield(KE_PARK(:,i),'raan');
end
surf(SMA,INC,raan)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('RAAN (deg)')
title('RAAN vs. SMA and INC')
hold off




end