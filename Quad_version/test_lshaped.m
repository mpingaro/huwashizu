%% TEST MESH
clear all; close all; clc;

ndx = 10;
ndy = 10;
length = 2;
heigth = 2;

[coordinates,element,npoint,nelem] = lshaped(length,heigth,ndx,ndy);

% PLOT SOLUTION 
for i =1:nelem
    x(1,[1,2,3,4]) = coordinates(element(i,[1,2,3,4]),1);
    y(1,[1,2,3,4]) = coordinates(element(i,[1,2,3,4]),2);
    % Undeformed mesh
    hold on
    surf(x,y,zeros(4))
end 
hold off    
axis equal
view(0,90)