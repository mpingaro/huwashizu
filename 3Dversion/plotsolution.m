% by Marco Pingaro & Paolo Venini

function plotsolution(coordinates,element,defo,strain,stress)


%% PLOT SOLUTION
% Undeformed mesh

figure
% Plot mesh
tetramesh(element, coordinates, ones(size(element,1),1),'FaceAlpha',1);
camorbit(20,0);
title('Undeformed Mesh','fontsize',14);
axis equal
% Deformed mesh
figure
% Plot mesh
tetramesh(element, defo, ones(size(element,1),1),'FaceAlpha',1);
camorbit(20,0);
title('Deformed Mesh','fontsize',14);
axis equal

%% STRAIN
strain = reshape(strain',6,[]);
strain_xx = strain(1,:);
strain_yy = strain(2,:);
strain_zz = strain(3,:);
strain_xy = strain(4,:);
strain_xz = strain(5,:);
strain_yz = strain(6,:);

% Media degli strain
st_xx = zeros(size(element,1),1);
st_yy = zeros(size(element,1),1);
st_zz = zeros(size(element,1),1);
st_xy = zeros(size(element,1),1);
st_xz = zeros(size(element,1),1);
st_yz = zeros(size(element,1),1);
for i=1:size(element,1)
    st_xx(i) = sum(strain_xx(element(i,:)))/4;
    st_yy(i) = sum(strain_yy(element(i,:)))/4;
    st_zz(i) = sum(strain_zz(element(i,:)))/4;
    st_xy(i) = sum(strain_xy(element(i,:)))/4;
    st_xz(i) = sum(strain_xz(element(i,:)))/4;
    st_yz(i) = sum(strain_yz(element(i,:)))/4;
end

figure
tetramesh(element, defo,st_xx);
camorbit(20,0);
title('STRAIN XX','fontsize',14);
axis equal

figure
tetramesh(element, defo,st_yy);
camorbit(20,0);
title('STRAIN YY','fontsize',14);
axis equal

figure
tetramesh(element, defo,st_zz);
camorbit(20,0);
title('STRAIN ZZ','fontsize',14);
axis equal

figure
tetramesh(element, defo,st_xy);
camorbit(20,0);
title('STRAIN XY','fontsize',14);
axis equal

figure
tetramesh(element, defo,st_xz);
camorbit(20,0);
title('STRAIN XZ','fontsize',14);
axis equal

figure
tetramesh(element, defo,st_yz);
camorbit(20,0);
title('STRAIN YZ','fontsize',14);
axis equal

%% STRESS
stress = reshape(stress',6,[]);
stress_xx = stress(1,:);
stress_yy = stress(2,:);
stress_zz = stress(3,:);
stress_xy = stress(4,:);
stress_xz = stress(5,:);
stress_yz = stress(6,:);

% Media degli strain
ss_xx = zeros(size(element,1),1);
ss_yy = zeros(size(element,1),1);
ss_zz = zeros(size(element,1),1);
ss_xy = zeros(size(element,1),1);
ss_xz = zeros(size(element,1),1);
ss_yz = zeros(size(element,1),1);
for i=1:size(element,1)
    ss_xx(i) = sum(stress_xx(element(i,:)))/4;
    ss_yy(i) = sum(stress_yy(element(i,:)))/4;
    ss_zz(i) = sum(stress_zz(element(i,:)))/4;
    ss_xy(i) = sum(stress_xy(element(i,:)))/4;
    ss_xz(i) = sum(stress_xz(element(i,:)))/4;
    ss_yz(i) = sum(stress_yz(element(i,:)))/4;
end

figure
tetramesh(element, defo,ss_xx);
camorbit(20,0);
title('STRESS XX','fontsize',14);
axis equal

figure
tetramesh(element, defo,ss_yy);
camorbit(20,0);
title('STRESS YY','fontsize',14);
axis equal

figure
tetramesh(element, defo,ss_zz);
camorbit(20,0);
title('STRESS ZZ','fontsize',14);
axis equal

figure
tetramesh(element, defo,ss_xy);
camorbit(20,0);
title('STRESS XY','fontsize',14);
axis equal

figure
tetramesh(element, defo,ss_xz);
camorbit(20,0);
title('STRESS XZ','fontsize',14);
axis equal

figure
tetramesh(element, defo,ss_yz);
camorbit(20,0);
title('STRESS YZ','fontsize',14);
axis equal

end