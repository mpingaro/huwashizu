% ----------------------------------------------------------------------- %
% ------ A FINITE ELEMENT METHOD FOR A THREE-FIELD FORMULATION OF ------- %  
% ------    LINEAR ELASTICITY BASED ON BIORTHOGONAL SYSTEMS       ------- %
% ----------------------------------------------------------------------- %
% ----------------   By (Marco Pingaro, Paolo Venini)   ----------------- %
% ----------------              3D VERSION              ----------------- %
%                                                                         %
% - Marco Pingaro    : Phd Student of University of Pavia                 %
%               mail : marco.pingaro@iusspavia.it                         %
% - Paolo Venini     : Professor of University of Pavia                   %
%               mail : paolo.venini@unipv.it                              %
% ----------------------------------------------------------------------- %
% SPACE APPROXIMATION:
% ------- Displacement : P1 + B  C^0               
% ------- Strain       : P1      C^0              
% ------- Stress       : P1      C^-1
% -- B is the space of boubble function
% ------------------------------------------------------------------------%
clear all; close all; clc;

global coordinates
%% INPUT
length = 5;
width  = 1;
thick  = 1;
ndx = 10;                                     % partition in x direction
ndy = 4;                                      % partition in y direction
ndz = 4;                                      % partition in z direction
young = 1.e3;                                 % young modulus
poisson = 0.3;                                % poisson modulus

% Neumann boudary conditions (edges) (not implemented)
bcn = [];                                     % index of edges 
fn(1,:) = [0, 0, 0];                          % Traction edge 1
fn(2,:) = [0, 0, 0];                          % Traction edge 2
fn(3,:) = [0, 0, 0];                          % Traction edge 3
fn(4,:) = [0, 0, 0];                          % Traction edge 4
fn(5,:) = [0, 0, 0];                          % Traction edge 5
fn(6,:) = [0, 0, 0];                          % Traction edge 6
% Neumann boudary conditions (nodes)
bct = [6,8];                              % index of edges 
ft(1,:) = [0, 0, 0];                          % Traction nodes 1 (0,0,0)
ft(2,:) = [0, 0, 0];                          % Traction nodes 2 (lx,0,0)
ft(3,:) = [0, 0, 0];                          % Traction nodes 3 (0,ly,0)
ft(4,:) = [0, 0, 0];                          % Traction nodes 4 (lx,ly,0)
ft(5,:) = [0, 0, 0];                          % Traction nodes 5 (0,0,lz)
ft(6,:) = [0, 0,-2];                          % Traction nodes 6 (lx,0,lz)
ft(7,:) = [0, 0, 0];                          % Traction nodes 7 (0,ly,lz)
ft(8,:) = [0, 0,-2];                          % Traction nodes 8 (lx,ly,lz)
% Dirichlet boudary conditions
bcd = 1 ;                                     % index of edges
ud(1,:) = [0, 0, 0];                          % Displacement edge 1 (x=0) 
ud(2,:) = [0, 0, 0];                          % Displacement edge 2 (x=lx)
ud(3,:) = [0, 0, 0];                          % Displacement edge 3 (z=0)
ud(4,:) = [0, 0, 0];                          % Displacement edge 4 (z=lz)
ud(5,:) = [0, 0, 0];                          % Displacement edge 5 (y=0)
ud(6,:) = [0, 0, 0];                          % Displacement edge 6 (y=ly)
%% Body load
g  = [0, 0, 0];                               % Body load

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
spost = solve_HuWashizu(KASSEM,F,length,width,thick,bcn,fn,bct,ft,bcd,ud);

%% POST PROCESSING
[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

%% PLOT SOLUTION
plotsolution(coordinates,element,defo,strain,stress);