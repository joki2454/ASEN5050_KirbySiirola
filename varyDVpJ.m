%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   VpJ1_jc
%   Dv_ramrange
%   SET
%   
% Outputs:
%   Dv_ramvec
%   a_park
%   i_park
%   
function [Dv_ramvec,a_park,i_park] = varyDVpJ(VpJ1_jc,Dv_ramrange,SET)
%% Produce ram vector Dvs and ram unit vector
uRamvec   = unitvec(VpJ1_jc); % km/s, ram unit vector
Dv_ramvec = linspace(Dv_ramrange(1),Dv_ramrange(2),SET.PRESENT.numSteps);

%% Solve for SMA and INC over the full DV range
for i = 1:length(Dv_ramvec)
  [a_park(i),i_park(i),TOF_JSOI,~,~,~,~] = transferSequence(Dv_ramvec(i).*uRamvec,SET);
end

end