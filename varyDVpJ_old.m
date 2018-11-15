%%  Author: Joshua Kirby
%  Created: 11/11/2018
% Modified: 11/11/2018
%
% Purpose: 
%
% Inputs:
%   VpJ1_jc
%   Dv_ramrange
%   Dv_rhatrange
%   SET
%   
% Outputs:
%   Dv_ramvec
%   Dv_rhatvec
%   a_park
%   i_park
%   
function [Dv_ramvec,Dv_rhatvec,DVpS,a_park,i_park] = varyDVpJ(RpJ_jc,VpJ1_jc,Dv_ramrange,Dv_rhatrange,SET)
%% Produce ram vector Dvs and ram unit vector
uRamvec   = unitvec(VpJ1_jc); % ram unit vector
Dv_ramvec = linspace(Dv_ramrange(1),Dv_ramrange(2),SET.PRESENT.numSteps);

%% Produce rhat vector Dvs and rhat unit vector
uRhatvec  = unitvec(RpJ_jc); % rhat unit vector
Dv_rhatvec= linspace(Dv_rhatrange(1),Dv_rhatrange(2),SET.PRESENT.numSteps);

%% Solve for SMA and INC over the full DV range
for i = 1:length(Dv_ramvec)
  for j = 1:length(Dv_rhatvec)
    [a_park(i,j),i_park(i,j),~,~,~,~,DVpS(:,i,j),~] = transferSequence(Dv_ramvec(i).*uRamvec + Dv_rhatvec(j).*uRhatvec,SET);
  end
end

end