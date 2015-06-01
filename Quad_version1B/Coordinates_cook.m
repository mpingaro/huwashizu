% by Marco Pingaro & Paolo Venini

function [coordinates,nnod]=Coordinates_cook(nodes,ndx,ndy,dl1,dl2)

%% Coorinates
% theta
nnod = (ndx+1)*(ndy+1);
theta = zeros(ndy+1,1);
for i = 0:length(theta)-1
    theta(i+1,1) = atand((nodes(2,2)+i*dl1/ndy-i*dl2/ndy)/nodes(2,1));
end

%% Coordinates Matrix
% This matrix contains at rows the nodal points of physical element and at 
% columns the relative coordinates
coordinates = zeros(nnod,2);
for i = 1:ndy+1
    for j = 1:ndx+1
        coordinates(j+(i-1)*(ndx+1),1) = (j-1)*nodes(2,1)/ndx;
        coordinates(j+(i-1)*(ndx+1),2) = coordinates(j+(i-1)*(ndx+1),1)*tand(theta(i,1))...
        +(i-1)*nodes(4,2)/ndy;
    end
end

return
