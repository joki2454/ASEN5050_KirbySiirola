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
function [] = presentResults(Dv_ramvec,Dv_rhatvec,DVpS,a_park,i_park,SET)
%% Determine total DeltaV for each Dv_ramvec
for i = 1:length(Dv_ramvec)
  for j = 1:length(Dv_rhatvec)
    Dv_tot(i,j) = norm([Dv_ramvec(i) Dv_rhatvec(j)]) + norm(DVpS(:,i,j)); % km/s
  end
end

figure
hold on
grid on
grid minor
surf(a_park,i_park,Dv_tot)
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('Total DeltaV (km/s)')
hold off




























end