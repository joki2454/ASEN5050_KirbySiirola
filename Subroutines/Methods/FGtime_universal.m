%%  Author: Joshua Kirby
%  Created: 11/12/2018
% Modified: 11/12/2018
%
% Purpose:
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
v0 = norm(v0); % km/s
ksi = v0^2/2-mu/r0;
alpha = -v0^2/mu+2/r0;

%% Calculate Orbital Elements
[a,e,~,~,~,~,~] = inertial2keplerian(R0,V0,mu);

%% Determine Alpha
alpha = mu/a; % km^2/s^2

%% Determine inital guess for chi
if alpha > 0.000001 % elliptical or circular
  chi0 = sqrt(mu)*TOF*alpha;
elseif abs(alpha) < 0.000001 % parabolic
  H = cross(R0,V0); % km^2/s
  p = norm(H)^2/mu; % km
  s = 1/2*acot(3*sqrt(mu/p^3)*TOF); % rad
  w = atan(tan(s)^(1/3)); % rad
  chi0 = sqrt(p)*2*cot(2*w);
elseif alpha < -0.000001 % hyperbolic
  A = 1/alpha;
  chi0 = sign(TOF)*sqrt(-A)*log(-2*mu*alpha*TOF/(dot(R0,V0)+sign(TOF)*sqrt(-mu*A)*(1-r0*alpha));
else
  error('Invalid value for alpha')
end

%% Solve for 
    
    
    
    
    
    
    
    
    

end