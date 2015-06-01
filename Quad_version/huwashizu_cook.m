% By Marco Pingaro & Paolo Venini

function sp = huwashizu_cook(ndx,ndy)

young = 250;                                  % young modulus
poisson = 0.4999;                             % poisson modulus
nodes   = [0, 0; 48, 44; 48, 60; 0, 44] ;
dl1     = nodes(3,2)-nodes(2,2) ;
dl2     = nodes(4,2) ;
% Neumann boudary conditions (edges)
bcn = 4;                                      % index of edges 
fn(1,:) = [0, 0];                             % Traction edge 1
fn(2,:) = [0, 0];                             % Traction edge 2
fn(3,:) = [0, 0];                             % Traction edge 3
fn(4,:) = [0, 100];                           % Traction edge 4
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
g       = [0, 0];                             % Body load

%% GEOMETRY
% Mesh 
[coordinates,nnod]=Coordinates_cook(nodes,ndx,ndy,dl1,dl2);
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
%[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

%% PLOT SOLUTION
%plotsol(coordinates,defo,strain,stress,ndx,ndy);

sp = max(spost);

return