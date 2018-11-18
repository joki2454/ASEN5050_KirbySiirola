%%  Author: Joshua Kirby
%  Created: 11/16/2018
% Modified: 11/16/2018
%
% Purpose: 
%
% Inputs:
%   R - position
%   r - radius
%   c - (1,3) color vector
%
function [h] = plotBody3D(R,r,c)
%% Produce 50x50 sphere of radius r at position R
[x,y,z] = sphere(50);
x = r*x + R(1);
y = r*y + R(2);
z = r*z + R(3);

%% Plot the body
for i = 1:3
  C(:,:,i) = ones(size(x))*c(i);
end
h = surf(x,y,z,C,'linestyle','none');

end