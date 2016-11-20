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
title('Deformation','fontsize', 24);
axis equal
view(0,90)
%print -color -depsc defo_beam.eps;

% STRAIN
strain = reshape(strain',3,[]);
strain_xx = strain(1,:);
strain_xy = strain(2,:);
strain_yy = strain(3,:);

strain_xx = reshape(strain_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xx)
title('Strain d_{xx}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_xx_beam.eps;

strain_xy = reshape(strain_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xy)
title('Strain d_{xy}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_xy_beam.eps;

strain_yy = reshape(strain_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_yy)
title('Srain d_{yy}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_yy_beam.eps;

% STRESS
stress = reshape(stress',3,[]);
stress_xx = stress(1,:);
stress_xy = stress(2,:);
stress_yy = stress(3,:);

stress_xx = reshape(stress_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xx)
title('Stress \sigma_{xx}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_xx_beam.eps;

stress_xy = reshape(stress_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xy)
title('Stress \sigma_{xy}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_xy_beam.eps;

stress_yy = reshape(stress_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_yy)
title('Stress \sigma_{yy}','fontsize', 24, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_yy_beam.eps;


% Save figure
print( figure(2),'-deps','-color','defo_beam.eps' );
print( figure(3),'-deps','-color','strain_xx_beam.eps' );
print( figure(4),'-deps','-color','strain_xy_beam.eps' );
print( figure(5),'-deps','-color','strain_yy_beam.eps' );
print( figure(6),'-deps','-color','stress_xx_beam.eps' );
print( figure(7),'-deps','-color','stress_xy_beam.eps' );
print( figure(8),'-deps','-color','stress_yy_beam.eps' );

end
