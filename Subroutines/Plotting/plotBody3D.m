%%  Author: Joshua Kirby
%  Created: 11/16/2018
% Modified: 11/16/2018
%
% Purpose: 
%
% Inputs:
%   R - position
%   r - radius
%
function [] = plotBody3D(R,r)
%% Produce 50x50 sphere of radius r at position R
[x,y,z] = sphere(50);
x = r*x + R(1);
y = r*y + R(2);
z = r*z + R(3);

%% Plot the body
surf(x,y,z)

end