%%  Author: Joshua Kirby
%  Created: 10/28/18
% Modified: 11/15/18
%
% Purpose: This function determines the initial DeltaV and final DeltaV
% required to follow a Lambert's arc transfer in a two-body problem from 
% an initial position and velocity to a final position and velocity given
% some time of flight.  A selection can be made between using a transfer
% with transfer angle less than 180 deg or greater than 180 deg.  Default
% is less than 180 deg.
%
% NOTES: -Transfer assumes transfer angle is less than 360 deg.
%        -All vectors must be given in consistent inertial coordinates. 
%
% Inputs (all units must be consistent): 
%   R1        - 3D initial position vector
%   V1        - 3D initial velocity vector
%   R2        - 3D final position vector
%   V2        - 3D final velocity vector
%   TOF       - Desired time of flight from initial to final position
%   mu        - gravitational parameter of system (or of central body if
%               spacecraft mass is significantly less than mass of central body)
%
% Varargin:
%   Name             Value               Description
%   ------------------------------------------------------------------------
%   transferAngle    'nu<180'            default, change in true anomaly of
%                                        transfer is less than 180 deg
%   transferAngle    'nu>180'            change in true anomaly of
%                                        transfer is greater than 180 deg
%   ------------------------------------------------------------------------
%
% Outputs (units consistent with inputs):
%   DV1       - 3D initial DeltaV vector
%   DV2       - 3D final DeltaV vector
%
function [DV1,DV2] = lambert_givenTOF(R1,V1,R2,V2,TOF,mu,varargin)
%% Operational parameters
tol = 1000*eps; % solution tolerance for semi-major axis

%% Determine Initial Transfer Properties
r1 = norm(R1);
r2 = norm(R2);

%% Process varargin
% defaults
nusign = 1; % sign of transfer angle (positive yields nu<180, negative yields nu>180)

% overwrite defaults with varargin
if ~isempty(varargin)
  for i = 1:length(varargin)
    if strcmpi(varargin{i},'transferAngle')
      if strcmpi(varargin{i+1},'nu<180')
        nusign = 1;
      elseif strcmpi(varargin{i+1},'nu>180')
        nusign = -1;
      else
        error('Invalid option for ''transferAngle'' varargin, use either ''nu<180'' or ''nu>180''');
      end % transferAngle internal
    end % transferAngle detect
  end % varargin for
end % varargin exists

%% Determine Transfer Angle
switch nusign
  case 1
    Dnu = nusign*abs(acosd(dot(R1,R2)/(r1*r2))); % deg
  case -1
    Dnu = nusign*abs(acosd(dot(R1,R2)/(r1*r2)))+360; % deg
  otherwise
    error('Invalid nusign! Fix this code, it''s broke.')
end

%% Determine geometric quanitities
c = sqrt(r1^2+r2^2-2*r1*r2*cosd(Dnu)); % chord between initial and final position
s = 0.5*(r1+r2+c);                     % semi-perimeter of transfer triangle

%% Determine if transfer is elliptic or hyperbolic
% calculate parabolic TOF
switch nusign
  case 1 % Dnu < 180 deg
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)-(s-c)^(3/2));
  case -1 % Dnu > 180 deg
    TOFp = 1/3*sqrt(2/mu)*(s^(3/2)+(s-c)^(3/2));
  otherwise
    error('Invalid nusign! Fix this code, it''s broke.')
end

% Error catching
if TOF == TOFp
  error('Transfer is perfectly parabolic, so this function doesn''t know what to do...')
end

% Determination
isElliptic = TOF>TOFp; % true if elliptic, false if hyperbolic

%% Find minimum energy TOF
am     = s/2;                   % minimum energy semi-major axis
nm     = sqrt(mu/am^3);         % minimum energy mean motion, rad/<time>^2
alpham = pi;                    % minimum energy alpha, rad
betam0 = 2*asin(sqrt((s-c)/s)); % minimum energy principal beta, rad
betam  = nusign*betam0;         % minimum energy beta, rad
TOFm   = 1/nm*((alpham-betam)-(sin(alpham)-sin(betam)));
alphasign = sign(TOFm-TOF);

%% Determine alpha and beta
if isElliptic
  alpha0 = @(a) 2.*asin(sqrt(s./(2.*a))); % rad
  beta0  = @(a) 2.*asin(sqrt((s-c)./(2.*a))); % rad
  beta   = @(a) nusign.*beta0(a); % rad
  switch alphasign
    case 1
      alpha = @(a) alpha0(a); % rad
    case -1
      alpha = @(a) 2.*pi-alpha0(a); % rad
  end
else % hyperbolic
  alpha = @(a) 2.*asinh(sqrt(s./(2.*abs(a)))); % rad
  beta  = @(a) 2.*asinh(sqrt((s-c)./(2.*abs(a)))); % rad
end

%% Formulate TOF equation
if isElliptic
  TOFeqn = @(a) sqrt(a.^3./mu).*(alpha(a)-beta(a)-(sin(alpha(a))-sin(beta(a)))) - TOF;
else % hyperbolic
  switch nusign
    case 1 % Dnu < 180
      TOFeqn = @(a) sqrt(abs(a).^3./mu).*(sinh(alpha(a))-alpha(a)-(sinh(beta(a))-beta(a))) - TOF;
    case -1 % Dnu > 180
      TOFeqn = @(a) sqrt(abs(a).^3./mu).*(sinh(alpha(a))-alpha(a)+(sinh(beta(a))-beta(a))) - TOF;
  end
end

%% Solve for semi-major axis
options = optimoptions('fsolve','maxiter',500,'tolfun',tol,'Display','none');
a = fsolve(TOFeqn,1.2*am,options);

% Error catching
if isnan(a)
  error('fsolve could not solve for semi-major axis');
end

if ~isElliptic % hyperbolic
  a = -abs(a);
end

%% Check validity of result
if abs(TOFeqn(a)) > 100*tol
  error('Semi-major axis solution does not produce desired TOF to within realistic tolerance');
end

%% Calculate transfer orbit parameters e, nu
% eccentricity
if isElliptic
  e = sqrt(1-4/c^2*(s-r1)*(s-r2)*sin((alpha(a)+beta(a))/2)^2);
else % hyperbolic
  e = sqrt(1+4/c^2*(s-r1)*(s-r2)*sinh((alpha(a)+beta(a))/2)^2);
end

% true anomalies
nutmp = acosd(1/e*(a*(1-e^2)/r1-1)); % deg
nu1s = [nutmp, 360-nutmp]; % deg
nutmp = acosd(1/e*(a*(1-e^2)/r2-1)); % deg
nu2s = [nutmp, 360-nutmp]; % deg

% Find combination of true anomalies whose difference is Dnu
for i = 1:length(nu1s)
  for j = 1:length(nu2s)
    if abs(nu2s(j)-nu1s(i)-Dnu)<10000*eps
      nu2 = nu2s(j); % deg
      nu1 = nu1s(i); % deg
    end
  end
end

%% Calculate initial and final velocity along transfer
p = a*(1-e^2);
% using f and g functions
f = 1-r2/p*(1-cosd(Dnu));
g = r1*r2/(sqrt(mu*p))*sind(Dnu);
fdot = sqrt(mu/p)*tand(Dnu/2)*((1-cosd(Dnu))/p-1/r2-1/r1);
gdot = 1 - r1/p*(1-cosd(Dnu));

% Find velocities
Vt1 = (R2-f*R1)/g;
Vt2 = fdot*R1+gdot*Vt1;

%% Calculate DeltaVs
DV1 = Vt1-V1;
DV2 = V2-Vt2;


end