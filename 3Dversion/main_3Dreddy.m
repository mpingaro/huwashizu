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
% ------- Displacement : vertor (1,2) = P1 + B  C^0               
% ------- Strain       : tensor (2,2) = P1      C^0              
% ------- Stress       : tensor (2,2) = P1      C^-1
% -- B is the space of boubble function
% ------------------------------------------------------------------------%
clear all; close all; clc;

%% INPUT
length = 5;
width  = 1;
thick  = 1;
ndx = 2;                                      % partition in x direction
ndy = 2;                                      % partition in y direction
ndz = 2;                                      % partition in z direction
young = 1.e3;                                 % young modulus
poisson = 0.3;                                % poisson modulus

% Neumann boudary conditions (edges)
bcn = [];                                     % index of edges 
fn(1,:) = [0, 0, 0];                          % Traction edge 1
fn(2,:) = [0, 0, 0];                          % Traction edge 2
fn(3,:) = [0, 0, 0];                          % Traction edge 3
fn(4,:) = [0, 0, 0];                          % Traction edge 4
fn(5,:) = [0, 0, 0];                          % Traction edge 5
fn(6,:) = [0, 0, 0];                          % Traction edge 6
% Neumann boudary conditions (nodes)
bct = [];                                     % index of edges 
ft(1,:) = [0, 0, 0];                          % Traction nodes 1
ft(2,:) = [0, 0, 0];                          % Traction nodes 2
ft(3,:) = [0, 0, 0];                          % Traction nodes 3
ft(4,:) = [0, 0, 0];                          % Traction nodes 4
ft(5,:) = [0, 0, 0];                          % Traction nodes 5
ft(6,:) = [0, 0, 0];                          % Traction nodes 6
ft(7,:) = [0, 0, 0];                          % Traction nodes 7
ft(8,:) = [0, 0, 0];                          % Traction nodes 8
% Dirichlet boudary conditions
bcd = [];                                     % index of edges
ud(1,:) = [0, 0, 0];                          % Displacement edge 1 
ud(2,:) = [0, 0, 0];                          % Displacement edge 2
ud(3,:) = [0, 0, 0];                          % Displacement edge 3
ud(4,:) = [0, 0, 0];                          % Displacement edge 4
ud(5,:) = [0, 0, 0];                          % Displacement edge 5
ud(6,:) = [0, 0, 0];                          % Displacement edge 6
%% Body load
g  = [0, 0, -1];                               % Body load

%% GEOMETRY
% Mesh 
[element,coordinates] = regular_tetrahedral_mesh(length,width,thick,ndx,ndy,ndz);
%
nnod  = size(coordinates,1);
nelem = size(element,1);
ngdls = 6*nnod;
%
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);  % Corrispondence Matrix displacement 
mc2=CorrispoMC2(element,nelem);             % Corrispondence Matrix strain, stress
%
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
plotsolution(coordinates,element,defo,strain,stress);
