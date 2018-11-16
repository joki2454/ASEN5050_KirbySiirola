%%  Author: Joshua Kirby
%  Created: 11/03/2018
% Modified: 11/15/2018
%
% Purpose: This function carries out an interplanetary transfer by beginning
% with Cassini's initial conditions just after it's Earth flyby (20 Aug
% 1999 00:00:00.000 UTC), entering a Jupiter flyby, performing a maneuver
% at perijove, exiting a Jupiter flyby, entering a Saturn insertion, and
% performing an antivelocity impulsive maneuver at perisaturnium to enter 
% a circular parking orbit.  This is carried out using the method of
% patched conics such that:
%   - from the initial time (denoted T0) to the time of
%     entering Jupiter's sphere-of-influence (denoted JSOI1) the spacecraft is
%     in a sun-spacecraft two-body problem (2BP)
%   - from JSOI1 to perijove (denoted PJ) the spacecraft is in a Jupiter-
%     spacecraft 2BP
%   - from PJ to exiting Jupiter's sphere-of-influence (denoted JSOI2) the 
%     spacecraft is in a Jupiter-spacecraft 2BP
%   - from JSOI2 to entering Saturn's sphere-of-influence (denoted SSOI) the 
%     spacecraft is in a sun-spacecraft 2BP
%   - from SSOI to perisaturnium (denoted PS) the spacecraft is in a
%     Saturn-spacecraft 2BP
% All solutions are analytical (not numerically integrated) since the
% spacecraft is always in a 2BP.  Within a planet's sphere-of-influence (SOI),
% calculations are performed based on the conic equation and Kepler's
% equation.  When traveling in a sun-spacecraft 2BP to reach Saturn's SOI
% the time of flight (TOF) to reach the SOI is determined via a targeter which varies the
% TOF until the final position is within a tolerance of Saturn's SOI.  
% The products of this function are parking orbit keplerian elements and all 
% relevant TOFs, positions, and velocities.
%
% NOTES:
%   - Vectors are expressed in a J2000 frame with various origins, which are
%     specified
%   - Parking orbit keplerian elements are expressed relative to a
%     Saturn-centered J2000 frame
%   - The TOF from T0 to SSOI1 was obtained via numerical integration in
%     STK and is simply used and verified here.  This was done to reduce
%     computational time for this function.  It is possible because the
%     only variation in the transfer sequence (DeltaV at perijove) is 
%     performed AFTER the flight from T0 to JSOI1.
%
% Assumptions: The planet's SOI's are assumed to be constant and to have a
% value equal to the exact value of the SOI for that planet at 00:00:00.000 UTC on
% Cassini's true flyby date (Jupiter: 30 Dec 2000, Saturn: 01 Jul 2004).
%
% Inputs: 
%   DVpJ           - 3D vector of DeltaV applied at Perijove, km/s
%   SET            - struct of settings, initial state, and options
%
% Outputs:
%   KE_PARK        - struct of keplerian elements of parking orbit:
%                    SMA (km), ECC (), INC (deg), AOP (deg), RAAN (deg)
%   TOF            - struct of various times of flight for transfer:
%                    JSOI       - time spent in Jupiter SOI, s
%                    JSOI2_SSOI - time spent between Jupiter SOI and
%                                     Saturn SOI, s
%                    SSOI       - time spent in Saturn SOI, s
%   T0             - struct of positions and velocities at initial time:
%                    R_hc - position in heliocentric    J2000 frame, km
%                    V_hc - velocity in heliocentric    J2000 frame, km/s
%   JSOI1          - struct of positions and velocities upon entering
%                       Jupiter SOI:
%                    R_jc - position in Jupiter-centric J2000 frame, km
%                    V_jc - velocity in Jupiter-centric J2000 frame, km/s
%                    R_hc - position in heliocentric    J2000 frame, km
%                    V_hc - velocity in heliocentric    J2000 frame, km/s
%   PJ             - struct of positions and velocities at perijove:
%                    R_jc  - position in Jupiter-centric J2000 frame, km
%                    V1_jc - velocity in Jupiter-centric J2000 frame before maneuver, km/s
%                    V2_jc - velocity in Jupiter-centric J2000 frame after  maneuver, km/s
%   JSOI2          - struct of positions and velocities upon exiting
%                       Jupiter SOI:
%                    R_jc - position in Jupiter-centric J2000 frame, km
%                    V_jc - velocity in Jupiter-centric J2000 frame, km/s
%                    R_hc - position in heliocentric    J2000 frame, km
%                    V_hc - velocity in heliocentric    J2000 frame, km/s           
%   SSOI           - struct of positions and velocities upon entering
%                       Saturn SOI:
%                    R_sc - position in Saturn-centric J2000 frame, km
%                    V_sc - velocity in Saturn-centric J2000 frame, km/s
%                    R_hc - position in heliocentric   J2000 frame, km
%                    V_hc - velocity in heliocentric   J2000 frame, km/s
%   PS             - struct of positions and velocities at perisaturnium:
%                    R_sc  - position in Saturn-centric   J2000 frame, km
%                    V1_sc - velocity in Saturn-centric   J2000 frame before maneuver, km/s
%                    V2_sc - velocity in Saturn-centric   J2000 frame after  maneuver, km/s
%                    DV    - maneuver at perisaturnium in J2000 frame, km/s
%   optim_badness  - use this?
%
function [KE_PARK,TOF,T0,JSOI1,PJ,JSOI2,SSOI,PS,optim_badness]...
              = transferSequence(DVpJ,SET)
%% Data extraction (heliocentric J2000 initial state)
T0.R_hc = SET.CASS.R0; % km
T0.V_hc = SET.CASS.V0; % km/s

%% Use FG Functions to propagate to Jupiter SOI in heliocentric J2000 2BP
[JSOI1.R_hc,JSOI1.V_hc] = FGtime_universal(T0.R_hc,T0.V_hc,SET.CONST.muSun,SET.CASS.TOF_0_JSOI1);

%% Calculate Jupiter state at t_JSOI1
jstate_JSOI1 = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1,'J2000','NONE','Sun');

%% Calculate spacecraft position and velocity in Jupiter-centered J2000 frame
JSOI1.R_jc = JSOI1.R_hc-jstate_JSOI1.state(1:3); % km
JSOI1.V_jc = JSOI1.V_hc-jstate_JSOI1.state(4:6); % km/s

%% Confirm that JSOI1 position is at Jupiter SOI
if abs(norm(JSOI1.R_jc)-SET.CONST.JSOI)/SET.CONST.JSOI*100 > 0.01
  error('Calculated position is not within 0.01% of the sphere of influence of Jupiter')
end

%% Calculate perijove position and velocity in J2000 frame using Jupiter-centered 2BP
% orbital elements and angular momentum magnitude at JSOI1
[a,e,i,omega,Omega,nuJSOI1,h] = inertial2keplerian(JSOI1.R_jc,JSOI1.V_jc,SET.CONST.muJ);
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
PJ.R_jc = C*RpJ_rth; % km,   position at perijove
PJ.V1_jc = C*VpJ_rth; % km/s, velocity at perijove before maneuver

% Time past periapsis at JSOI1
E = TAtoEA(e,nuJSOI1);
t_tpJSOI1 = sqrt(a^3/SET.CONST.muJ)*(E-e*sin(E));

%% Perform impulsive maneuver at perijove
PJ.V2_jc = PJ.V1_jc+DVpJ; % km/s, velocity at perijove after maneuver

%% Calculate JSOI2 position and velocity in J2000 frame using Jupiter-centered 2BP
% orbital elements and angular momentum magnitude at JSOI2
[a,e,i,omega,Omega,~,h] = inertial2keplerian(PJ.R_jc,PJ.V2_jc,SET.CONST.muJ);

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
JSOI2.R_jc = C*RJSOI2_rth; % km
JSOI2.V_jc = C*VJSOI2_rth; % km/s

% Time past periapsis at JSOI2
E = TAtoEA(e,nu);
t_tpJSOI2 = sqrt(a^3/SET.CONST.muJ)*(E-e*sin(E));

%% Determine JSOI2 position and velocity in J2000 heliocentric frame
% time spent in Jupiter's sphere of influence
TOF.JSOI = t_tpJSOI2-t_tpJSOI1; % s

% Jupiter's state at JSOI2 in heliocentric J2000 frame
jstate_JSOI2 = mice_spkezr('Jupiter Barycenter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1+TOF.JSOI,'J2000','NONE','Sun');

% J2000 heliocentric position and velocity at JSOI2
JSOI2.R_hc = JSOI2.R_jc + jstate_JSOI2.state(1:3); % km
JSOI2.V_hc = JSOI2.V_jc + jstate_JSOI2.state(4:6); % km/s

%% Use targeter to find TOF from JSOI2 to SSOI
% Epoch at JSOI2
etJSOI2 = cspice_str2et(SET.CASS.startDate) + SET.CASS.TOF_0_JSOI1 + TOF.JSOI; % s

% Formulate targeter equation
SSOI_targeter_eqn = @(TOF) SSOI_targeter(JSOI2.R_hc,JSOI2.V_hc,SET.CONST.muSun,TOF,etJSOI2,SET);

% Solve for TOF, using STK (DeltaV=0) TOF from JSOI2 to SSOI as starting point
options = optimset('TolFun',1e-04,'display','none');
[TOF.JSOI2_SSOI,fval,exitflag,~] = fminsearch(SSOI_targeter_eqn,SET.CASS.TOF_JSOI2_SSOI,options); % s

%% Error Checking
if ismember(exitflag,-1:0)
  error('fminsearch produced an error')
end
  
%% Define result badness
if fval <= sqrt(options.TolFun)
  optim_badness = 0; % Cassini hit SSOI
elseif fval > sqrt(options.TolFun)
  optim_badness = 1; % Cassini never hit SSOI
else
  error('Fminsearch has produced an invalid fval, scary')
end

%% Determine SSOI position and velocity in J2000 heliocentric frame
[SSOI.R_hc,SSOI.V_hc] = FGtime_universal(JSOI2.R_hc,JSOI2.V_hc,SET.CONST.muSun,TOF.JSOI2_SSOI);

%% Determine SSOI position and velocity in J2000 Saturn-centered frame
% Saturn state at SSOI
sstate_SSOI = mice_spkezr('Saturn Barycenter',etJSOI2+TOF.JSOI2_SSOI,'J2000','NONE','Sun');

% J2000 saturn-centered position and velocity at SSOI
SSOI.R_sc = SSOI.R_hc - sstate_SSOI.state(1:3); % km
SSOI.V_sc = SSOI.V_hc - sstate_SSOI.state(4:6); % km/s

%% Confirm that SSOI position is at Saturn SOI
% THIS ERROR IS NOW CAPTURED BY RESULT BADNESS (fsolve_badness)

if abs(norm(SSOI.R_sc)-SET.CONST.SSOI)/SET.CONST.SSOI*100 > 0.01
  % error('Calculated position is not within 0.01% of the sphere of influence of Saturn')
end

%% Calculate perisaturnium position and velocity in J2000 frame using Saturn-centered 2BP
% orbital elements and angular momentum magnitude at SSOI
[a,e,i,omega,Omega,nuSSOI,h] = inertial2keplerian(SSOI.R_sc,SSOI.V_sc,SET.CONST.muS);
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
PS.R_sc = C*RpS_rth; % km,   position at perisaturnium
PS.V1_sc = C*VpS_rth; % km/s, velocity at perisaturnium before maneuver

% Time past periapsis at SSOI
E = TAtoEA(e,nuSSOI);
t_tpSSOI = sqrt(a^3/SET.CONST.muS)*(E-e*sin(E));

%% Calculate semi-major axis of circular parking orbit
a_park = norm(PS.R_sc); % km

%% Calculate DeltaV at perisaturnium required to enter circular parking orbit
% Determine magnitude of required DeltaV at perisaturnium (in ram direction)
DvpS = sqrt(SET.CONST.muS/a_park) - norm(PS.V1_sc);

% Determine DeltaV vector at perisaturnium in Saturn-centered J2000 frame
PS.DV = DvpS.*unitvec(PS.V1_sc); % km/s

% Time past periapsis at perisaturnium
t_tppS = 0; % s

% Time spent in Saturn's SOI
TOF.SSOI = t_tppS - t_tpSSOI; % s

%% Determine orbital elements of parking orbit
% Parking orbit velocity at perisaturnium in Saturn-centered frame
PS.V2_sc = PS.V1_sc + PS.DV; % km/s

% Calculate orbital elements
[KE_PARK.a,KE_PARK.e,KE_PARK.i,KE_PARK.aop,KE_PARK.raan,~,~] = inertial2keplerian(PS.R_sc,PS.V2_sc,SET.CONST.muS);

end