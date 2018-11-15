%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/12/2018
%
% Purpose:  
%
% Inputs:
%   
%   
% Outputs:
%   
function [] = Copy_of_presentResults(DVpJ,SET)
%% 
sma = linspace(SET.RANGES.sma(1),SET.RANGES.sma(2),SET.PRESENT.numSteps); % sma
inc = linspace(SET.RANGES.inc(1),SET.RANGES.inc(2),SET.PRESENT.numSteps); % deg
for i = 1:length(sma)
  for j = 1:length(inc)
    Dv_tot(i,j) = norm(squeeze(DVpJ(i,j,:))); % km/s
  end
end

[SMA,INC] = meshgrid(sma,inc);

figure
hold on
grid on
grid minor
surf(SMA,INC,Dv_tot)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('Total DeltaV (km/s)')
hold off




























end