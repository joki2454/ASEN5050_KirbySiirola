%%  Author: Joshua Kirby
%  Created: 10/28/2018
% Modified: 10/30/2018
%
% Purpose: This function will convert an inertial position and velocity
% vector of a spacecraft into the keplerian orbital elements of that craft
% in the same inertial system.
%
% TODO: Implement a tolerance rather than checking for exact zeros
% throughout.
%
% Inputs (all must be in consistent units):
%   R - position vector
%   V - velocity vector
%   mu - gravitational parameter of central body
%   varargin - 'nupos' true anomaly will range from [0 360] (default)
%              'nuneg' true anomaly will range from [-180 180]
%
% Outputs (same units as input, all angles are degrees)
%   a - semi-major axis
%   e - eccentricity (magnitude)
%   i - inclination
%   omega - argument of perigee [variable range]
%   Omega - right ascension of the ascending node
%   nu - true anomaly
%   h - specific angular momentum (magnitude)
function [a,e,i,omega,Omega,nu,h] = inertial2keplerian(R,V,mu,varargin)
%% Ensure that inputs are column vectors
if size(R) == [3 1]
  % all good
elseif size(R) == [1 3]
  R = R';
else
  error('R must be a 3-element vector')
end
if size(V) == [3 1]
  % all good
elseif size(V) == [1 3]
  V = V';
else
  error('V must be a 3-element vector')
end
%% Process varargin
% Defaults
nusign = 1;

% Update defaults with varargin
if ~isempty(varargin)
  for i = 1:length(varargin)
    if strcmpi(varargin{i},'nupos')
      nusign = 1;
    elseif strcmpi(varargin{i},'nuneg')
      nusign = -1;
    else
      error('Invalid varargin, allowed values are ''nupos'' or ''nuneg''');
    end
  end
end

%% Helpful parameters
X = [1 0 0]';
Y = [0 1 0]';
Z = [0 0 1]';
r = norm(R);
v = norm(V);
H = cross(R,V); % angular momentum vector
h = norm(H);
epsilon = v^2/2-mu/r; % specific energy
N = cross(Z,H); % line of nodes
n = norm(N);
E = cross(V,H)/mu-R/r; % eccentricity vector

%% Orbital elements (with sign checks)
i     = acosd(H(3)/h);
a     = -mu/(2*epsilon);
e     = sqrt(1+2*h^2*epsilon/mu^2);
Omega = sign(dot(N,Y))*abs(acosd(dot(N,X)/n));
omega = sign(dot(E,Z))*abs(acosd(dot(N,E)/(n*e)));
nu    = sign(dot(R,V))*abs(acosd(dot(R,E)/(r*e)));

%% Correct range of true anomaly
switch nusign
  case 1
    nu = nu + 360;
    while nu > 360
      nu = nu - 360;
    end
  case -1
    nu = nu; % all is good
  otherwise
    error('Invalid nusign.  Fix this code, it''s broke.');
end














