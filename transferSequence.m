%%  Author: Joshua Kirby
%  Created: 11/03/2018
% Modified: 11/12/2018
%
% Purpose: 
%
% Inputs: 
%   DVpJ           - DeltaV applied at Perijove, km/s
%   SET            - struct of settings, initial state, and options
%
% Outputs:
%   a_park         - semi-major axis of final circular parking orbit, km
%   i_park         - inclination of final circular parking orbit, deg
%   TOF_JSOI       - time spent in Jupiter's SOI, s
%   TOF_JSOI2_SSOI - time spent between Saturn's and Jupiter's SOI, s
%   TOF_SSOI       - time spent in between Saturn's SOI and perisaturnium, s
%   RpJ_jc         - position at perijove in Jupiter-centered J2000 frame,
%                    km
%   VpJ1_jc        - velocity at perijove before a maneuver in a
%                    Jupiter-centered J2000 frame, km/s
%   DVpS           - DeltaV maneuver at perisaturnium in J2000 frame to 
%                    put spacecraft into a circular parking orbit, km/s
%   optim_badness  - 0 - all is well
%                    1 - Solution not found in max iterations of fsolve or an error occurred, this 
%                        could mean that 1) max iterations is
%                        too small, 2) there is some problem with the fsolve
%                        setup, or 3) (purpose of this variable) satellite does 
%                        not pass through SSOI
%
function [a_park,i_park,TOF_JSOI,TOF_JSOI2_SSOI,TOF_SSOI,...
    RpJ_jc,VpJ1_jc,DVpS,optim_badness] = transferSequence(DVpJ,SET)
%% Data extraction (heliocentric J2000 initial state)
R0 = SET.CASS.R0; % km
V0 = SET.CASS.V0; % km/s

%% Use FG Functions to propagate to Jupiter SOI in heliocentric J2000 2BP
[RJSOI1_hc,VJSOI1_hc] = FGtime_universal(R0,V0,SET.CONST.muSun,SET.CASS.TOF_0_JSOI1);

%% Calculate Jupiter state at t_JSOI1
jstate_JSOI1 = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1,'J2000','NONE','Sun');

%% Calculate spacecraft position and velocity in Jupiter-centered J2000 frame
RJSOI1_jc = RJSOI1_hc-jstate_JSOI1.state(1:3); % km
VJSOI1_jc = VJSOI1_hc-jstate_JSOI1.state(4:6); % km/s

%% Confirm that JSOI1 position is at Jupiter SOI
if abs(norm(RJSOI1_jc)-SET.CONST.JSOI)/SET.CONST.JSOI*100 > 0.01
  error('Calculated position is not within 0.01% of the sphere of influence of Jupiter')
end

%% Calculate perijove position and velocity in J2000 frame using Jupiter-centered 2BP
% orbital elements and angular momentum magnitude at JSOI1
[a,e,i,omega,Omega,nuJSOI1,h] = inertial2keplerian(RJSOI1_jc,VJSOI1_jc,SET.CONST.muJ);
nu = 0; % deg, true anomaly is 0 deg at perijove

% Distance from origin at periapsis
r = a*(1-e^2)/(1+e*cosd(nu));

% Rotating frame velocity components
vr = SET.CONST.muJ/h*e*sind(nu); % km/s
vt = SET.CONST.muJ/h*(1+e*cosd(nu)); % km/s

% Formulate vectors in rotating frame
RpJ_rth = [r 0 0]'; % km
VpJ_rth = [vr vt 0]'; % km/s

% Rotation from rotating to J2000 frame
C = Crth2xyz(i,omega,Omega,nu);

% Position and velocity vectors in jupiter-centered J2000 frame
RpJ_jc = C*RpJ_rth; % km,   position at perijove
VpJ1_jc = C*VpJ_rth; % km/s, velocity at perijove before maneuver

% Time past periapsis at JSOI1
E = TAtoEA(e,nuJSOI1);
t_tpJSOI1 = sqrt(a^3/SET.CONST.muJ)*(E-e*sin(E));

%% Perform impulsive maneuver at perijove
VpJ2_jc = VpJ1_jc+DVpJ; % km/s, velocity at perijove after maneuver

%% Calculate JSOI2 position and velocity in J2000 frame using Jupiter-centered 2BP
% orbital elements and angular momentum magnitude at JSOI2
[a,e,i,omega,Omega,~,h] = inertial2keplerian(RpJ_jc,VpJ2_jc,SET.CONST.muJ);

% Determine true anomaly at JSOI2 from conic equation,
%   sign must be positive because this is the outgoing branch of the hyperbola
nu   = abs(acosd(1/e*((a*(1-e^2))/SET.CONST.JSOI-1))); % deg

% distance from origin at JSOI2 is JSOI

% Rotating frame velocity components
vr = SET.CONST.muJ/h*e*sind(nu);
vt = SET.CONST.muJ/h*(1+e*cosd(nu));

% Formulate vectors in rotating frame
RJSOI2_rth = [SET.CONST.JSOI 0 0]'; % km
VJSOI2_rth = [vr vt 0]'; % km/s

% Rotation from rotating to J2000 frame
C = Crth2xyz(i,omega,Omega,nu);

% Position and velocity vectors in Jupiter-centered J2000 frame
RJSOI2_jc = C*RJSOI2_rth; % km
VJSOI2_jc = C*VJSOI2_rth; % km/s

% Time past periapsis at JSOI2
E = TAtoEA(e,nu);
t_tpJSOI2 = sqrt(a^3/SET.CONST.muJ)*(E-e*sin(E));

%% Determine JSOI2 position and velocity in J2000 heliocentric frame
% time spent in Jupiter's sphere of influence
TOF_JSOI = t_tpJSOI2-t_tpJSOI1; % s

% Jupiter's state at JSOI2 in heliocentric J2000 frame
jstate_JSOI2 = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1+TOF_JSOI,'J2000','NONE','Sun');

% J2000 heliocentric position and velocity at JSOI2
RJSOI2_hc = RJSOI2_jc + jstate_JSOI2.state(1:3); % km
VJSOI2_hc = VJSOI2_jc + jstate_JSOI2.state(4:6); % km/s

%% Use targeter to find TOF from JSOI2 to SSOI
% Epoch at JSOI2
etJSOI2 = cspice_str2et(SET.CASS.startDate) + SET.CASS.TOF_0_JSOI1 + TOF_JSOI; % s

% Formulate targeter equation
SSOI_targeter_eqn = @(TOF) SSOI_targeter(RJSOI2_hc,VJSOI2_hc,SET.CONST.muSun,TOF,etJSOI2,SET);

% Solve for TOF, using STK (DeltaV=0) TOF from JSOI2 to SSOI as starting point
options = optimset('TolFun',1e-04,'display','none');
[TOF_JSOI2_SSOI,fval,exitflag,~] = fminsearch(SSOI_targeter_eqn,SET.CASS.TOF_JSOI2_SSOI,options); % s

%% Error Checking
if ismember(exitflag,-1:0)
  error('fminsearch produced an error')
end
  
%% Define result badness
if fval < sqrt(options.TolFun)
  optim_badness = 0; % Cassini hit SSOI
elseif fval > sqrt(options.TolFun)
  optim_badness = 1; % Cassini never hit SSOI
else
  error('Fminsearch has produced an invalid fval, scary')
end

%% Determine SSOI position and velocity in J2000 heliocentric frame
[RSSOI_hc,VSSOI_hc] = FGtime_universal(RJSOI2_hc,VJSOI2_hc,SET.CONST.muSun,TOF_JSOI2_SSOI);

%% Determine SSOI position and velocity in J2000 Saturn-centered frame
% Saturn state at SSOI
sstate_SSOI = mice_spkezr('Saturn Barycenter',etJSOI2+TOF_JSOI2_SSOI,'J2000','NONE','Sun');

% J2000 saturn-centered position and velocity at SSOI
RSSOI_sc = RSSOI_hc - sstate_SSOI.state(1:3); % km
VSSOI_sc = VSSOI_hc - sstate_SSOI.state(4:6); % km/s

%% Confirm that SSOI position is at Saturn SOI
% THIS ERROR IS NOW CAPTURED BY RESULT BADNESS (fsolve_badness)

if abs(norm(RSSOI_sc)-SET.CONST.SSOI)/SET.CONST.SSOI*100 > 0.01
  % error('Calculated position is not within 0.01% of the sphere of influence of Saturn')
end

%% Calculate perisaturnium position and velocity in J2000 frame using Saturn-centered 2BP
% orbital elements and angular momentum magnitude at SSOI
[a,e,i,omega,Omega,nuSSOI,h] = inertial2keplerian(RSSOI_sc,VSSOI_sc,SET.CONST.muS);
nu = 0; % deg, true anomaly is 0 deg at perisaturnium

% Distance from origin at periapsis
r = a*(1-e^2)/(1+e*cosd(nu)); % km

% Rotating frame velocity components
vr = SET.CONST.muS/h*e*sind(nu); % km/s
vt = SET.CONST.muS/h*(1+e*cosd(nu)); % km/s

% Formulate vectors in rotating frame
RpS_rth = [r 0 0]'; % km
VpS_rth = [vr vt 0]'; % km/s

% Rotation from rotating to J2000 frame
C = Crth2xyz(i,omega,Omega,nu);

% Position and velocity vectors in Saturn-centered J2000 frame
RpS_sc = C*RpS_rth; % km,   position at perisaturnium
VpS1_sc = C*VpS_rth; % km/s, velocity at perisaturnium before maneuver

% Time past periapsis at SSOI
E = TAtoEA(e,nuSSOI);
t_tpSSOI = sqrt(a^3/SET.CONST.muS)*(E-e*sin(E));

%% Calculate semi-major axis of circular parking orbit
a_park = norm(RpS_sc); % km

%% Calculate DeltaV at perisaturnium required to enter circular parking orbit
% Determine magnitude of required DeltaV at perisaturnium (in ram direction)
DvpS = sqrt(SET.CONST.muS/a_park) - norm(VpS1_sc);

% Determine DeltaV vector at perisaturnium in Saturn-centered J2000 frame
DVpS = DvpS.*unitvec(VpS1_sc); % km/s

% Time past periapsis at perisaturnium
t_tppS = 0; % s

% Time spent in Saturn's SOI
TOF_SSOI = t_tppS - t_tpSSOI; % s

%% Determine inclination of circular parking orbit
% Parking orbit velocity at perisaturnium in Saturn-centered frame
VpS2_sc = VpS1_sc + DVpS; % km/s

% Angular momentum
h = cross(RpS_sc,VpS2_sc); % km^2/s

% Inclination
i_park = acosd(h(3)/norm(h)); % deg

end