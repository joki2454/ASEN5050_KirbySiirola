%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   
%   
% Outputs:
%   
function [] = presentResults(Dv_ramvec,a_park,i_park)
%% Temporary

figure
hold on
grid on
grid minor
plot3(a_park,i_park,Dv_ramvec,'o-')
xlabel('SMA (km)')
ylabel('INC (deg)')
zlabel('Total DeltaV (km/s)')
hold off




























end