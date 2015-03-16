% by Marco Pingaro & Paolo Venini

%function plotsolution_goal(coordinates,element,defo,strain,stress,X)
function plotsolution_goal(coordinates,element,X)

global nelem 

%[x1,x2] = size(X);
%XNEW = zeros(1,3*x1);
%for i = 1:nelem-1
%    XNEW(1,3*i+1) = X(i);
%    XNEW(1,3*i+2) = X(i);
%    XNEW(1,3*i+3) = X(i);
%end
%XNEW = 1-XNEW;
figure
colormap(gray)
ShowContDisp_goal(element,coordinates,X)
title('OPTIMAL TOPOLOGY','fontsize',14);
axis equal
view(0,90)

return

end