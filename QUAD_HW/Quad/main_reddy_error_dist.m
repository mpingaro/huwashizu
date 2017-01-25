% ----------------------------------------------------------------------- %
% ------ A FINITE ELEMENT METHOD FOR A THREE-FIELD FORMULATION OF ------- %  
% ------    LINEAR ELASTICITY BASED ON BIORTHOGONAL SYSTEMS       ------- %
% ----------------------------------------------------------------------- %
% ----------------   By (Marco Pingaro, Paolo Venini)   ----------------- %
%                                                                         %
% - Marco Pingaro    : Phd Student of University of Pavia                 %
%               mail : marco.pingaro@iusspavia.it                         %
% - Paolo Venini     : Professor of University of Pavia                   %
%               mail : paolo.venini@unipv.it                              %
% ----------------------------------------------------------------------- %
% SPACE APPROXIMATION:
% ------- Displacement : vertor (1,2) = Q1 + B   C^0               
% ------- Strain       : tensor (2,2) = Q1       C^0              
% ------- Stress       : tensor (2,2) = Q1       C^-1
% -- B is the single boubble functions
% ------------------------------------------------------------------------%
clear; close all; clc;

%% INPUT
length = 1;                                   % length 
height = 1;                                   % heigth
ndx = 1;                                      % partition in x direction
ndy = 1;                                      % partition in y direction
% Lamé constants
nu = 0.4999;
mu = 1;
%
lambda = 2*mu*nu/(1-2*nu);

% Neumann boudary conditions (edges)
bcn = [];                                     % index of edges 
fn(1,:) = [0, 0];                             % Traction edge 1
fn(2,:) = [0, 0];                             % Traction edge 2
fn(3,:) = [0, 0];                             % Traction edge 3
fn(4,:) = [0, 0];                             % Traction edge 4
% Neumann boudary conditions (vertex)
bct = [];                                     % index of edges 
ft(1,:) = [0, 0];                             % Traction vertex 1
ft(2,:) = [0, 0];                             % Traction vertex 2
ft(3,:) = [0, 0];                             % Traction vertex 3
ft(4,:) = [0, 0];                             % Traction vertex 4
% Dirichlet boudary conditions
bcd = [1,2,3,4];                              % index of edges
ud(1,:) = [0, 0];                             % Displacement edge 1 
ud(2,:) = [0, 0];                             % Displacement edge 2
ud(3,:) = [0, 0];                             % Displacement edge 3
ud(4,:) = [0, 0];                             % Displacement edge 4
% Alpha coefficient

cf = [1,2,3]; % coefficient
nl = [4, 8, 16, 32, 64];
% Cycle on different coeffients
for k=1:3

alpha = cf(k)*mu;

name =   sprintf('elastic_error_disp_u_l2_1B_type1_%dmu_dist.txt',cf(k));
name_1 = sprintf('elastic_error_l2_strain_xx_1B_type1_%dmu_dist.txt',cf(k));
name_2 = sprintf('elastic_error_l2_strain_yy_1B_type1_%dmu_dist.txt',cf(k));
name_3 = sprintf('elastic_error_l2_strain_xy_1B_type1_%dmu_dist.txt',cf(k));

f = fopen( name, 'w' );
f1 = fopen(name_1, 'w');
f2 = fopen(name_2, 'w');
f3 = fopen(name_3, 'w');

fprintf(f, 'elements v.s. error in L2 norm\n'); 
fprintf(f1,'elements v.s. error in L2 norm strain xx\n');
fprintf(f2,'elements v.s. error in L2 norm strain yy\n');
fprintf(f3,'elements v.s. error in L2 norm strain xy\n');

% Cycle on element
for i=1:size(nl,2)

ndx = nl(i);
ndy = nl(i);

%% GEOMETRY
dx = length/ndx;
dy = height/ndy;
% Mesh 
%[coordinates,nnod]=Coordinates(ndx,ndy,dx,dy);
coordinates = beam_distorted(length, height, ndx, ndy);
nnod = size(coordinates,1);
[element,nelem]=Element(ndx,ndy);
%
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);  % Corrispondence Matrix displacement 
mc2=CorrispoMC2(element,nelem);             % Corrispondence Matrix strain, stress

ngdls = 3*nnod;
%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM,F,D,W,B,M,K] = assembly_error(coordinates,element,mc,mc2,lambda,alpha,mu,nelem,ngdlu,ngdls);

%% SOLVE
spost = solve_HuWashizu(KASSEM,F,ndx,ndy,bcn,fn,bct,ft,bcd,ud);

%% POST PROCESSING
[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

% Compute error in norm L2 displacement
er_u = error_l2_norm(spost, coordinates, lambda);    

% Compute error in norm L2 strain
[er_exx,er_eyy, er_exy] = error_strain_l2_norm(strain, coordinates, lambda);

% Print results
fprintf(f,  '%6.0f \t %6.5e \n', nelem, er_u);
fprintf(f1, '%6.0f \t %6.5e \n', nelem, er_exx);
fprintf(f2, '%6.0f \t %6.5e \n', nelem, er_eyy);
fprintf(f3, '%6.0f \t %6.5e \n', nelem, er_exy);

end
fclose(f);
fclose(f1);
fclose(f2);
fclose(f3);
end
