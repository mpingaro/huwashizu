% by Marco Pingaro & Paolo Venini

function plotsol(coordinates,defo,ndx,ndy)

% PLOT SOLUTION 
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;

defo_x = reshape(defo(:,1),ndx+1,ndy+1) ;
defo_y = reshape(defo(:,2),ndx+1,ndy+1) ;

% Undeformed mesh
figure, surf(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Deformed mesh
figure, surf(defo_x,defo_y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

end
