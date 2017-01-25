% by Marco Pingaro & Paolo Venini

function plotsol(coordinates,defo,strain,stress,ndx,ndy)

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

% STRAIN
strain = reshape(strain',3,[]);
strain_xx = strain(1,:);
strain_xy = strain(2,:);
strain_yy = strain(3,:);

strain_xx = reshape(strain_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xx)
title('STRAIN XX','fontsize',14);
axis equal
view(0,90)

strain_xy = reshape(strain_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xy)
title('STRAIN XY','fontsize',14);
axis equal
view(0,90)

strain_yy = reshape(strain_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_yy)
title('STRAIN YY','fontsize',14);
axis equal
view(0,90)

% STRESS
stress = reshape(stress',3,[]);
stress_xx = stress(1,:);
stress_xy = stress(2,:);
stress_yy = stress(3,:);

stress_xx = reshape(stress_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xx)
title('STRESS XX','fontsize',14);
axis equal
view(0,90)

stress_xy = reshape(stress_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xy)
title('STRESS XY','fontsize',14);
axis equal
view(0,90)

stress_yy = reshape(stress_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_yy)
title('STRESS YY','fontsize',14);
axis equal
view(0,90)

end
