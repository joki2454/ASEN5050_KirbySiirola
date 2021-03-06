%%  Author: Joshua Kirby
%  Created: 11/16/2018
% Modified: 11/17/2018
%
% Purpose: 
%
% Inputs:
%   DVpJ - 
%   SET  - 
%
% Outputs:
%   t       - (:) vector of times, s
%   R_jc    - (3, )
%   R_sc    - (3, )
%   R_hc    - (3,:) array of 3D spacecraft positions in heliocentric J2000 frame, km
%   P_EPHEM - struct of specific planetary ephemeris useful for plotting
function [t,R_hc,R_jc,R_sc,P_EPHEM] = transferStates(TOF,T0,JSOI1,PJ,JSOI2,SSOI,PS,KE_PARK,SET)
%% Generate time vector
periodPark = 2*pi*sqrt(KE_PARK.a^3/SET.CONST.muS); % period of parking orbit around saturn
t0     = 0;                                        % s
tJSOI1 = t0 + SET.CASS.TOF_0_JSOI1;                % s
tpJ    = tJSOI1 + TOF.JSOI1_PJ;                    % s
tJSOI2 = tpJ + TOF.PJ_JSOI2;                       % s
tSSOI  = tJSOI2 + TOF.JSOI2_SSOI;                  % s
tpS    = tSSOI + TOF.SSOI;                         % s
tend   = tpS + 1*periodPark;                       % s
t      = [linspace(t0,tJSOI1,    round(SET.PRESENT.timeSteps/10))...
          linspace(tJSOI1,tJSOI2,round(SET.PRESENT.timeSteps/10))...
          linspace(tJSOI2,tSSOI, round(SET.PRESENT.timeSteps/10))...
          linspace(tSSOI,tend,   round(SET.PRESENT.timeSteps*7/10))];  % s

%% Generate planetary ephemeris
jstate_pJ     = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate) + tpJ,'J2000','NONE','Sun');
P_EPHEM.RJ_pJ = jstate_pJ.state(1:3); % km
P_EPHEM.VJ_pJ = jstate_pJ.state(4:6); % km/s

sstate_pS     = mice_spkezr('Saturn Barycenter',cspice_str2et(SET.CASS.startDate) + tpS,'J2000','NONE','Sun');
P_EPHEM.RS_pS = sstate_pS.state(1:3); % km
P_EPHEM.VS_pS = sstate_pS.state(4:6); % km/s

%% Produce positions
cJ = 1; % Jupiter-centric counter
cS = 1; % Saturn- centric counter
for i = 1:length(t)
  % Produce planetary positions
  jstate = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate)+t(i),'J2000','NONE','Sun');
  sstate = mice_spkezr('Saturn Barycenter',cspice_str2et(SET.CASS.startDate)+t(i),'J2000','NONE','Sun');
  P_EPHEM.RJ_hc(:,i) = jstate.state(1:3); % km
  P_EPHEM.RS_hc(:,i) = sstate.state(1:3); % km
  if (t(i) >= t0) && ( t(i) < tJSOI1)
    % From T0 to JSOI1
    Dt = t(i) - t0; % s
    R_hc(:,i) = FGtime_universal(T0.R_hc,T0.V_hc,SET.CONST.muSun,Dt); % km
    
  elseif ( t(i) >= tJSOI1) && ( t(i) < tpJ)
    % From JSOI1 to PJ
    Dt = t(i) - tJSOI1; % s
    R_jc(:,cJ) = FGtime_universal(JSOI1.R_jc,JSOI1.V_jc,SET.CONST.muJ,Dt); % km
    R_hc(:,i)  = R_jc(:,cJ) + P_EPHEM.RJ_hc(:,i); % km
    cJ = cJ + 1;
    
  elseif ( t(i) >= tpJ) && ( t(i) < tJSOI2)
    % From PJ to JSOI2
    Dt = t(i) - tpJ; % s
    R_jc(:,cJ) = FGtime_universal(PJ.R_jc,PJ.V2_jc,SET.CONST.muJ,Dt); % km
    R_hc(:,i)  = R_jc(:,cJ) + P_EPHEM.RJ_hc(:,i); % km
    cJ = cJ + 1;
    
  elseif ( t(i) >= tJSOI2) && ( t(i) < tSSOI)
    % While between JSOI and SSOI
    Dt = t(i) - tJSOI2; % s
    R_hc(:,i) = FGtime_universal(JSOI2.R_hc,JSOI2.V_hc,SET.CONST.muSun,Dt); % km
    
  elseif ( t(i) >= tSSOI) && ( t(i) < tpS)
    % While within SSOI
    Dt = t(i) - tSSOI; % s
    R_sc(:,cS) = FGtime_universal(SSOI.R_sc,SSOI.V_sc,SET.CONST.muS,Dt); % km
    R_hc(:,i)  = R_sc(:,cS) + P_EPHEM.RS_hc(:,i);  %km
    cS = cS + 1;
    
  elseif ( t(i) >= tpS ) && ( t(i) < tend)
    % Saturn parking orbit
    Dt = t(i) - tpS; % s
    R_sc(:,cS) = FGtime_universal(PS.R_sc,PS.V2_sc,SET.CONST.muS,Dt); % km
    R_hc(:,i)  = R_sc(:,cS) + P_EPHEM.RS_hc(:,i); % km
    cS = cS + 1;
    
  else
    % do nothing ?
  end
end












end