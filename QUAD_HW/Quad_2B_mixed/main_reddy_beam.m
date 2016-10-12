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
% ------- Displacement : vertor (1,2) = Q1 + B1 + B2  C^0               
% ------- Strain       : tensor (2,2) = Q1            C^0              
% ------- Stress       : tensor (2,2) = Q1            C^-1
% -- B1 and B2 (mixed) are the two boubble functions
% ------------------------------------------------------------------------%
clear all; close all; clc;

%% INPUT
length = 10;                                  % length 
height = 2;                                   % heigth
nx = [4,8,16,32,64,128];                      % partition in x direction
ny = [2,4, 8,16,32,64] ;                      % partition in y direction
young = 1500;                                 % young modulus
poisson = 0.4999;                             % poisson modulus
f = 300;                                      % max value of distributed load
cf = [1,2,3];


for k=1:numel(cf)
fname = sprintf('error_beam_u_l2_%dmu.txt',cf(k));
f = fopen( fname, 'w');
fprintf(f, 'element vs. error u in norm L2\n');

for i=1:numel(nx)

ndx = nx(i);
ndy = ny(i);
%% GEOMETRY
dx = length/ndx;
dy = height/ndy;
% Mesh 
[coordinates,nnod]=Coordinates(ndx,ndy,dx,dy);
[element,nelem]=Element(ndx,ndy);
%
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);  % Corrispondence Matrix displacement 
mc2=CorrispoMC2(element,nelem);             % Corrispondence Matrix strain, stress

ngdls = 3*nnod;
lambda = young*poisson/( (1+poisson)*(1-2*poisson) );
mu = young/(2*(1+poisson));
alpha = cf(k)*mu;

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM,D,W,B,M,K] = assembly_beam(coordinates,element,mc,mc2,lambda,alpha,mu,nelem,ngdlu,ngdls);

%% SOLVE
spost = solve_HuWashizu_beam(KASSEM,coordinates,height,f,ndx,ndy,ngdlu);

%% POST PROCESSING
[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

%% Computing L2 error
er_u = error_beam_l2_norm(spost, coordinates, height, young, poisson, f);
fprintf(f, '%6.0f \t %6.5e \n', nelem, er_u);

end
fclose(f);
%% PLOT SOLUTION
%plotsol(coordinates,defo,strain,stress,ndx,ndy)
end
