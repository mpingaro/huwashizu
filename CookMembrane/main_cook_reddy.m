% ----------------------------------------------------------------------- %
% ------ A FINITE ELEMENT METHOD FOR A THREE-FIELD FORMULATION OF ------- %  
% ------    LINEAR ELASTICITY BASED ON BIORTHOGONAL SYSTEMS       ------- %
% ----------------------------------------------------------------------- %
% ------           By (Marco Pingaro & Paolo Venini)              ------- %
%                                                                         %
% - Marco Pingaro : Phd Student of University of Pavia                    %
%                   Dipartimento di Ingegneria Civile e Architettura      %
%            mail : marco.pingaro@iusspavia.it                            %
% - Paolo Venini  : Professor of University of Pavia                      %
%                   Dipartimento di Ingegneria Civile e Architettura      %
%            mail : paolo.venini@unipv.it                                 %
%                                                                         %
% ----------------------------------------------------------------------- %
% SPACE APPROXIMATION:
% ------- Displacement : vertor (1,2) = P1 + B  C^0               
% ------- Strain       : tensor (2,2) = P1      C^0              
% ------- Stress       : tensor (2,2) = P1      C^-1
% -- B is the space of boubble function!
% ------------------------------------------------------------------------%
%%
clear; close all; clc;
%% INPUT DATA
young = 250;
poisson = 0.4999;
f_traction = 100;

[ndx,ndy] = inputcook();
NODES = [0 0; 48 44; 48 60; 0 44]; 
DL1 = NODES(3,2)-NODES(2,2); 
DL2 = NODES(4,2);

%% COORDINATES MATRIX
[coordinates] = coordcook(NODES,ndx,ndy,DL1,DL2);
%% ELEMENT MATRIX
[element] = elcook(ndx,ndy);
%
nelem = size(element,1);
nnod  = size(coordinates,1);
ngdls = 3*nnod;
ngdlu = 2*nnod+2*nelem;

%% CORRESPONDENCE MATRIX (element)
[mc] = mccook(element);            % Corrispondence Matrix displacement
mc2=CorrispoMC2(element,nelem);    % Corrispondence Matrix strain, stress

lambda = young*poisson/( (1+poisson)*(1-2*poisson) );
mu = young/(2*(1+poisson));
alpha = 1;

%% ASSEMBLY GLOBAL SYSTEM
[KASSEM,F,D,W,B,M,K] = assembly(coordinates,element,mc,mc2,lambda,alpha,mu,nelem,ngdlu,ngdls);

%% SOLVE
[spost,defo,strain,stress] = solve_HuWashizu(coordinates,ndx,ndy,f_traction,KASSEM,ngdlu,D,W,B,M,K,alpha);

%% PLOT SOLUTION
%plotsolution(coordinates,element,defo,strain,stress);
