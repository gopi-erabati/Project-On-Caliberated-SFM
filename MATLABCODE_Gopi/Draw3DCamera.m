function [] = Draw3DCamera(p,cameraName,color)
% Draws the camera position and axis rotation in the 3D scene.
% 
% Input:
%	- p is a matrix at the form [r c] where
%            c is the coordinates of the camera, and
%	         r is the rotation matrix.
%
%	- cameraName is a text displayed at left of the camera. It is optional.
%
%   - color is the color of z-axis (principal axis).
%

if (nargin < 3)
    color = 'r';
end

c = p(:,4);

% Plot a circle on the 3D coordinates of the camera
plot3(c(1),c(3),c(2),'o','MarkerSize',10,'Color',color,'Linewidth',1);
%axis equal;