% by Marco Pingaro & Paolo Venini

function plotsol_lshaped(coordinates,element,defo)

nelem = size(element,1);
% PLOT SOLUTION 

% Undeformed mesh
figure,
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

figure,
% Deformed mesh
for i =1:nelem
    dx(1,[1,2]) = defo(element(i,[1,2]),1);
    dx(2,[1,2]) = defo(element(i,[4,3]),1);
    dy(1,[1,2]) = defo(element(i,[1,2]),2);
    dy(2,[1,2]) = defo(element(i,[4,3]),2);
    % Undeformed mesh
    hold on
    surf(dx,dy,zeros(2))
end 
hold off    
axis equal
view(0,90)

% % STRAIN
% strain = reshape(strain',3,[]);
% strain_xx = strain(1,:);
% strain_xy = strain(2,:);
% strain_yy = strain(3,:);
% 
% strain_xx = reshape(strain_xx,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,strain_xx)
% title('STRAIN XX','fontsize',14);
% axis equal
% view(0,90)
% 
% strain_xy = reshape(strain_xy,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,strain_xy)
% title('STRAIN XY','fontsize',14);
% axis equal
% view(0,90)
% 
% strain_yy = reshape(strain_yy,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,strain_yy)
% title('STRAIN YY','fontsize',14);
% axis equal
% view(0,90)
% 
% % STRESS
% stress = reshape(stress',3,[]);
% stress_xx = stress(1,:);
% stress_xy = stress(2,:);
% stress_yy = stress(3,:);
% 
% stress_xx = reshape(stress_xx,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,stress_xx)
% title('STRESS XX','fontsize',14);
% axis equal
% view(0,90)
% 
% stress_xy = reshape(stress_xy,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,stress_xy)
% title('STRESS XY','fontsize',14);
% axis equal
% view(0,90)
% 
% stress_yy = reshape(stress_yy,ndx+1,ndy+1) ;
% figure
% surf(defo_x,defo_y,stress_yy)
% title('STRESS YY','fontsize',14);
% axis equal
% view(0,90)

end