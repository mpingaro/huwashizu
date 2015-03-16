% by Marco Pingaro & Paolo Venini

function plotsolution_andmises(coordinates,element,defo,strain,stress,x)

global qexp SY 

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
hold off
colorbar
view(0,90)

%
figure
ShowContDisp(element,defo,strain_xy)
title('STRAIN XY','fontsize',14);
axis equal
hold off
colorbar
view(0,90)
%
figure
ShowContDisp(element,defo,strain_yy)
title('STRAIN YY','fontsize',14);
axis equal
hold off
colorbar
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
hold off
colorbar
view(0,90)
%
figure
ShowContDisp(element,defo,stress_xy)
title('STRESS XY','fontsize',14);
axis equal
hold off
colorbar
view(0,90)
%
figure
ShowContDisp(element,defo,stress_yy)
title('STRESS YY','fontsize',14);
axis equal
hold off
colorbar
view(0,90)

%
smises2 = stress_xx.^2+stress_yy.^2-stress_xx.*stress_yy+3*stress_xy.^2;
figure
ShowContDisp(element,defo,smises2.^(1/2))
title('MISES STRESS','fontsize',14);
hold off
colorbar
axis equal
view(0,90)

%stress evaluation
% STRESS
stress = reshape(stress',3,[]);
sxx = stress(1,:)';
sxy = stress(2,:)';
syy = stress(3,:)';
%sforzi elementari
s_xx = 1/3*(sxx(element(:,1))+ sxx(element(:,2))+ sxx(element(:,3)));
s_xy = 1/3*(sxy(element(:,1))+ sxy(element(:,2))+ sxy(element(:,3)));
s_yy = 1/3*(syy(element(:,1))+ syy(element(:,2))+ syy(element(:,3)));  
smises2 = s_xx.^2+s_yy.^2-s_xx.*s_yy+3*s_xy.^2; 
fconst = smises2-SY^2*x.^(2*qexp);
fconst3 = ((smises2).^(1/2))./(SY*x.^(qexp));

imax = 0;
maxfc = -1.e7;

for ii=1:length(fconst3)
    if fconst3(ii) > maxfc
        maxfc = fconst3(ii);
        imax = ii;
    end
end

figure
%colormap(gray)
ShowContDisp_goal(element,coordinates,fconst)
title('Yield Check \sigma_{VM}^2 - \sigma_{Y}^2 x^{2q}','fontsize',14);
axis equal
hold off
colorbar
view(0,90)


end