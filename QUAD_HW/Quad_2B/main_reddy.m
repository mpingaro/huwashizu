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
% -- B1 and B2 are the two boubble functions
% ------------------------------------------------------------------------%
clear all; close all; clc;

%% INPUT
length = 4;                                   % length 
height = 2;                                   % heigth
ndx =  2;                                     % partition in x direction
ndy =  1;                                     % partition in y direction
young = 1.e3;                                 % young modulus
poisson = 0.2;                             % poisson modulus
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
bcd = 3;                                      % index of edges
ud(1,:) = [0, 0];                             % Displacement edge 1 
ud(2,:) = [0, 0];                             % Displacement edge 2
ud(3,:) = [0, 0];                             % Displacement edge 3
ud(4,:) = [0, 0];                             % Displacement edge 4
% Body load
g       = [0,-10];                            % Body load

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
alpha = 2*mu;

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM,F,D,W,B,M,K] = assembly(coordinates,element,mc,mc2,lambda,alpha,mu,g,nelem,ngdlu,ngdls);

%% SOLVE
spost = solve_HuWashizu(KASSEM,F,ndx,ndy,bcn,fn,bct,ft,bcd,ud);

%% POST PROCESSING
[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

%% PLOT SOLUTION
%plotsol(coordinates,defo,strain,stress,ndx,ndy)
