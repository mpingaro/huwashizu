% by Marco Pingaro & Paolo Venini

function plotsol_beam(coordinates,defo,strain,stress,ndx,ndy,height,young,poisson,ld)

% PLOT SOLUTION 
x = reshape(coordinates(:,1),ndx+1,ndy+1) ;
y = reshape(coordinates(:,2),ndx+1,ndy+1) ;

defo_x = reshape(defo(:,1),ndx+1,ndy+1) ;
defo_y = reshape(defo(:,2),ndx+1,ndy+1) ;

% Undeformed mesh
figure, surf(x,y,zeros(ndx+1,ndy+1))
axis equal
view(0,90)

% Deformed mesh (Computed)
figure, surf(defo_x,defo_y,zeros(ndx+1,ndy+1))
title('Deformation','fontsize', 20);
axis equal
view(0,90)
%print -color -depsc defo_beam.eps;

% Deformed mesh (Analytical)
an_def_x = 2*ld/(young*height)*(1-poisson^2).*coordinates(:,1).*(height/2-coordinates(:,2));
an_def_y = ld/(young*height).*( coordinates(:,1).^2 + poisson/(1-poisson).*( coordinates(:,2).^2 - height.*coordinates(:,2)));

an_defo_x = reshape(an_def_x,ndx+1,ndy+1) ;
an_defo_y = reshape(an_def_y,ndx+1,ndy+1) ;
figure, surf(an_defo_x,an_defo_y,zeros(ndx+1,ndy+1))
title('Deformation (Anlytical)','fontsize', 20);
axis equal
view(0,90)
%

% STRAIN
strain = reshape(strain',3,[]);
strain_xx = strain(1,:);
strain_xy = strain(2,:);
strain_yy = strain(3,:);

% Strain xx
strain_xx = reshape(strain_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xx)
title('Strain d_{xx}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_xx_beam.eps;

% Strain xx (Analytical)
an_strain_xx = 2*ld/(young*height)*(1-poisson^2).*(height/2-coordinates(:,2));
an_strain_xx = reshape(an_strain_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,an_strain_xx)
title('Strain d_{xx} (Anlytical)','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%

% Strain xy
strain_xy = reshape(strain_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_xy)
title('Strain d_{xy}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_xy_beam.eps;

% Strain yy
strain_yy = reshape(strain_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,strain_yy)
title('Srain d_{yy}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc strain_yy_beam.eps;

% Strain yy (Analytical)
an_strain_yy = ld/(young*height)*(poisson/(1-poisson)).*(2.*coordinates(:,2)-height);
an_strain_yy = reshape(an_strain_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,an_strain_yy)
title('Srain d_{yy} (Analytical)','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%
% STRESS
stress = reshape(stress',3,[]);
stress_xx = stress(1,:);
stress_xy = stress(2,:);
stress_yy = stress(3,:);

stress_xx = reshape(stress_xx,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xx)
title('Stress \sigma_{xx}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_xx_beam.eps;

stress_xy = reshape(stress_xy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_xy)
title('Stress \sigma_{xy}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_xy_beam.eps;

stress_yy = reshape(stress_yy,ndx+1,ndy+1) ;
figure
surf(defo_x,defo_y,stress_yy)
title('Stress \sigma_{yy}','fontsize', 20, 'interpreter', 'tex');
colorbar
axis equal
view(0,90)
%print -color -depsc stress_yy_beam.eps;


% Save figure
print( figure(2),'-deps','-color','defo_beam.eps' );
print( figure(3),'-deps','-color','defo_beam_analytical.eps' );
print( figure(4),'-deps','-color','strain_xx_beam.eps' );
print( figure(5),'-deps','-color','strain_xx_beam_analytical.eps' );
print( figure(6),'-deps','-color','strain_xy_beam.eps' );
print( figure(7),'-deps','-color','strain_yy_beam.eps' );
print( figure(8),'-deps','-color','strain_yy_beam_analytical.eps' );
print( figure(9),'-deps','-color','stress_xx_beam.eps' );
print( figure(10),'-deps','-color','stress_xy_beam.eps' );
print( figure(11),'-deps','-color','stress_yy_beam.eps' );

end
