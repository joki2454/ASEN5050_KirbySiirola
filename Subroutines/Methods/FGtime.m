%%  Author: Joshua Kirby
%  Created: 11/05/2018
% Modified: 11/05/2018
%
% Purpose:
%
% Req'd Routines:
%   -inertial2keplerian.m
%   -TAtoEA.m
%
% Inputs (shown in SI, but any consistent unit set will produce the same units):
%   R0   - 3D initial position vector, km
%   V0   - 3D initial velocity vector, km/s
%   mu   - gravitational parameter of central body, km^3/s^2
%   TOF  - time of flight, s
%
% Outputs:
%   R    - 3D final position vector, km
%   V    - 3D final velocity vector, km/s
%
function [R,V] = FGtime(R0,V0,mu,TOF)
%% Operational parameters
tol = 1000*eps; % fsolve tolerance for E, rad

%% Calculate DeltaE
% Obtain orbital elements
[a,e,~,~,~,nu0,~] = inertial2keplerian(R0,V0,mu);

% Calculate intial eccentric anomaly
E0 = TAtoEA(e,nu0); % rad

% Calculate mean motion
n = sqrt(mu/a^3); % rad/s

% Formulate kepler's equation
EAeqn = @(Evar) (Evar - e*sin(Evar)) - (E0 - e*sin(E0)) - n*TOF;

% Solve for final eccentric anomaly
options = optimoptions('fsolve','maxIter',500,'tolfun',tol,'Display','none');
E = fsolve(EAeqn,E0,options); % rad

% Calculate change in eccentric anomaly
DE = E-E0; % rad

%% Useful parameters
r0 = norm(R0); % km
v0 = norm(V0); % km/s

%% Formulate FG functions
f    = 1 - a/r0*(1-cos(DE));
g    = TOF - sqrt(a^3/mu)*(DE-sin(DE));

%% Solve for R
R = f*R0+g*V0; % km
r = norm(R); % km

%% Formulate FGdot functions
fdot = -sin(DE)*sqrt(mu*a)/(r0*r);
gdot = 1-a/r*(1-cos(DE));

%% Solve for V
V = fdot*R0+gdot*V0; % km/s

end