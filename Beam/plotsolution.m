% by Marco Pingaro & Paolo Venini

function plotsolution(coordinates,element,defo,strain,stress)


%% PLOT SOLUTION
% Undeformed mesh
figure
plotmesh(element,coordinates)
title('Undeformed Mesh','fontsize',14);
axis equal
% Deformed mesh
figure
plotmesh(element,defo)
title('Deformed Mesh','fontsize',14);
axis equal

% STRAIN
strain = reshape(strain',3,[]);
strain_xx = strain(1,:);
strain_xy = strain(2,:);
strain_yy = strain(3,:);

figure
ShowContDisp(element,defo,strain_xx)
title('STRAIN XX','fontsize',14);
axis equal
view(0,90)
%
figure
ShowContDisp(element,defo,strain_xy)
title('STRAIN XY','fontsize',14);
axis equal
view(0,90)
%
figure
ShowContDisp(element,defo,strain_yy)
title('STRAIN YY','fontsize',14);
axis equal
view(0,90)

% STRESS
stress = reshape(stress',3,[]);
stress_xx = stress(1,:);
stress_xy = stress(2,:);
stress_yy = stress(3,:);

figure
ShowContDisp(element,defo,stress_xx)
title('STRESS XX','fontsize',14);
axis equal
view(0,90)
%
figure
ShowContDisp(element,defo,stress_xy)
title('STRESS XY','fontsize',14);
axis equal
view(0,90)
%
figure
ShowContDisp(element,defo,stress_yy)
title('STRESS YY','fontsize',14);
axis equal
view(0,90)

end