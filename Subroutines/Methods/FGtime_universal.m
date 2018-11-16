%%  Author: Joshua Kirby
%  Created: 11/12/2018
% Modified: 11/15/2018
%
% Purpose:  Generated using Algorithm 8 in Vallado.  This function uses the
% universal variable formulation of f and g functions to find the final position
% and velocity of an object in a two-body problem given an initial position
% and velocity and a time of flight.
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
function [R,V] = FGtime_universal(R0,V0,mu,TOF)
%% Operational parameters
tol = 1000*eps; % fsolve tolerance for E, rad

%% Useful parameters
r0 = norm(R0); % km
v0 = norm(V0); % km/s
alpha = -v0^2/mu+2/r0;

%% Calculate SMA
a = 1/alpha; % km

%% Determine inital guess for chi
if alpha > eps % elliptical or circular
  chi = sqrt(mu)*TOF*alpha;
elseif abs(alpha) < eps % parabolic
  H = cross(R0,V0); % km^2/s
  p = norm(H)^2/mu; % km
  s = 1/2*acot(3*sqrt(mu/p^3)*TOF); % rad
  w = atan(tan(s)^(1/3)); % rad
  chi = sqrt(p)*2*cot(2*w);
elseif alpha < -eps % hyperbolic
  chi = sign(TOF)*sqrt(-a)*log(-2*mu*alpha*TOF/(dot(R0,V0)+sign(TOF)*sqrt(-mu*a)*(1-r0*alpha)));
else
  error('Invalid value for alpha')
end

%% Solve for psi and chi
solved = 0;
for i = 1:1000 % max iterations
  chilast = chi;
  psi = chilast^2*alpha;
  [c2,c3] = c2c3(psi);
  r = chilast^2*c2 + dot(R0,V0)/sqrt(mu)*chilast*(1-psi*c3)+r0*(1-psi*c2);
  chi = chilast + (sqrt(mu)*TOF-chi^3*c3-dot(R0,V0)/sqrt(mu)*chi^2*c2-r0*chi*(1-psi*c3))/r;
  if abs(chi-chilast) < 1e-06
    solved = 1;
    psi = chi^2*alpha;
    [c2,c3] = c2c3(psi);
    break
  end
end
    
if ~solved
  error('FGtime_universal failed to solve for chi to within tolerance in under max iterations')
end
    
%% Solve for f and g
f = 1 - chi^2/r0*c2;
g = TOF - chi^3/sqrt(mu)*c3;

%% Solve for R
R = f*R0 + g*V0; % km

%% Solve for fdot and gdot
r = norm(R); % km
fdot = sqrt(mu)/(r*r0)*chi*(psi*c3-1);
gdot = 1 - chi^2/r*c2;

%% Solve for V
V = fdot*R0 + gdot*V0; % km/s

end % FGtime_universal()



%% c2c3
% Generated using Algorithm 1 in Vallado
function [c2,c3] = c2c3(psi)
if psi > eps
  c2 = (1-cos(sqrt(psi)))/psi;
  c3 = (sqrt(psi)-sin(sqrt(psi)))/sqrt(psi^3);
else
  if psi < -eps
    c2 = (1-cosh(sqrt(-psi)))/psi;
    c3 = (sinh(sqrt(-psi))-sqrt(-psi))/sqrt((-psi)^3);
  else
    c2 = 1/2;
    c3 = 1/6;
  end
end
end