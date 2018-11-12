%%  Author: Joshua Kirby
%  Created: 11/03/2018
% Modified: 11/05/2018
%
% Purpose: 
%
% Inputs: 
%
% Outputs:
%
function [a_park,i_park,TOF_JSOI,TOF_JSOI2_SSOI] = transferSequence(DVpJ,SET)
%% Data extraction (heliocentric J2000 initial state)
R0 = SET.CASS.R0; % km
V0 = SET.CASS.V0; % km/s

%% Use FG Functions to propagate to Jupiter SOI in heliocentric J2000 2BP
[RJSOI1_hc,VJSOI1_hc] = FGtime(R0,V0,SET.CONST.muSun,SET.CASS.TOF_0_JSOI1);

%% Calculate Jupiter state after at t_JSOI1
jstate_JSOI1 = mice_spkezr('Jupiter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1,'J2000','NONE','Sun');

%% Calculate spacecraft position in Jupiter-centered J2000 frame
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
jstate_JSOI2 = mice_spkezr('Jupiter',cspice_str2et(SET.CASS.startDate)+SET.CASS.TOF_0_JSOI1+TOF_JSOI,'J2000','NONE','Sun');

% J2000 heliocentric position and velocity at JSOI2
RJSOI2_hc = RJSOI2_jc + jstate_JSOI2.state(1:3); % km
VJSOI2_hc = VJSOI2_jc + jstate_JSOI2.state(4:6); % km/s

%% Use targeter to find TOF from JSOI2 to SSOI
% Epoch at JSOI2
etJSOI2 = cspice_str2et(SET.CASS.startDate) + SET.CASS.TOF_0_JSOI1 + TOF_JSOI; % s

% Formulate targeter equation
SSOI_targeter_eqn = @(TOF) SSOI_targeter(RJSOI2_hc,VJSOI2_hc,SET.CONST.muSun,TOF,etJSOI2,SET);

% Solve for TOF, using STK (DeltaV=0) TOF from JSOI2 to SSOI as starting point
options = optimoptions('fsolve','maxiter',500,'tolfun',SET.TRGT.tol,'Display','none');
TOF_JSOI2_SSOI = fsolve(SSOI_targeter_eqn,SET.CASS.TOF_JSOI2_SSOI,options); % s

%% Determine SSOI position and velocity in J2000 heliocentric frame
[RSSOI_hc,VSSOI_hc] = FGtime(RJSOI2_hc,VJSOI2_hc,SET.CONST.muSun,TOF_JSOI2_SSOI);

%% Determine SSOI position and velocity in J2000 Saturn-centered frame
% Saturn state at SSOI
sstate_SSOI = mice_spkezr('Saturn',etJSOI2+TOF_JSOI2_SSOI,'J2000','NONE','Sun');

% J2000 saturn-centered position and velocity at SSOI
RSSOI_jc = RSSOI_hc - sstate.state(1:3); % km
VSSOI_jc = VSSOI_hc - sstate.state(4:6); % km/s

%% 






















end