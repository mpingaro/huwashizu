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
% ------- Displacement : vector (1,2) = P1 + B  C^0               
% ------- Strain       : tensor (2,2) = P1      C^0              
% ------- Stress       : tensor (2,2) = P1      C^-1
% -- B is the space of boubble function
% ------------------------------------------------------------------------%
clear all; close all; clc;

global YOUNG young youngmin pois coordinates element mc mc2 alpha g nelem ngdlu ngdls
global ndx ndy nnod
global AELE BELE WELE KELE MELE DELE FIXEDDOFS FORCEDDOFS FREEDOFS
global A B W M D FORIGINAL iut IOPTION IFILTER pexp HH strain stress spost 

%IOPTION = compliance come: [1 Lext; 2 C_epsi:epsi; 3 D_sigma:epsi
IOPTION = 2;
IFILTER = 0;
pexp = 3;

%% INPUT
length = 6;                                   % length 
height = 2;                                   % heigth
ndx = 6;                                     % partition in x direction
ndy = 2;                                     % partition in y direction
dx = length/ndx;
dy = height/ndy;
young = 10;                                    % young modulus virgin
youngmin = 1e-6;                              % minimum young modulus
poisson = 0.3;                                % poisson modulus
VOLFRAC = 0.35;
VMAX = length*height*VOLFRAC;

%% GEOMETRY
dx = length/ndx;
dy = height/ndy;
% Mesh 
[coordinates,nnod]=CoordinatesType2(ndx,ndy,dx,dy);
[element,nelem]=ElementType2(ndx,ndy);
%
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);  % Corrispondence Matrix displacement 
mc2=CorrispoMC2(element,nelem);             % Corrispondence Matrix strain, stress
iut = [1:4:nelem; 2:4:nelem; 3:4:nelem; 4:4:nelem]';

ngdls = 3*nnod;
lambda = young*poisson/( (1+poisson)*(1-2*poisson) );
mu = young/(2*(1+poisson));
alpha = 2*mu;
pois(1) = (1-poisson)/((1+poisson)*(1-2*poisson));
pois(2) = poisson/((1+poisson)*(1-2*poisson));
pois(3) = 1/(1+poisson);
YOUNG = young*ones(nelem,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GDL LATI
LX1 = 1:2:2*ndx+1; LY1 = LX1+1;
LX2 = 2*ndx+1:2*(ndx+1):2*(ndx+1)*(ndy+1)-1; LY2 = LX2+1;
LX3 = sort(2*(ndx+1)*(ndy+1)-1:-2:2*(ndx)*(ndy+1)-3); LY3 = LX3+1;
LX4 = 1:2*(ndx+1):2*(ndx)*(ndy+1)-3; LY4=LX4+1;

%GDL NODI
NX1 = min(LX1); NY1 = min(LY1);
NX2 = min(LX2); NY2 = min(LY2); 
NX3 = max(LX3); NY3 = max(LY3);
NX4 = max(LX4); NY4 = max(LY4);

FIXEDDOFS = [];
FIXEDDOFS(1,:) = [LX4 NY2]; FIXEDDOFS(2,:) = 0;

FORCEDDOFS = [];
FORCEDDOFS(1,:) = NY4; FORCEDDOFS(2,:) = -0.025;

% Body load
g = [0, 0];                            
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%element matrices with unitary Young
AELE = zeros(8,8,4); % A
BELE = zeros(9,8,4); % B
WELE = zeros(9,8,4); % W
KELE = zeros(9,9,4); % K
MELE = zeros(9,9,4); % M
DELE = zeros(9,9,4); % D
for k = 1:4
    P(1,[1 2]) = coordinates(element(k,1),[1 2]);
    P(1,[3 4]) = coordinates(element(k,2),[1 2]);
    P(1,[5 6]) = coordinates(element(k,3),[1 2]);
    [AELE(:,:,k),BELE(:,:,k),KELE(:,:,k),MELE(:,:,k),WELE(:,:,k),DELE(:,:,k)] = reddy_element(P,1);
end

%FILTER
if IFILTER == 1
    rmin = 1.*max(dx,dy);
    [bari,distanze] = distances(coordinates,element);
    HH = zeros(size(distanze));
    HH = max(0,rmin-distanze);
end

%assemblo una volta per tutte le matrici che non dipendono da Young
[A,B,W,M,D,FORIGINAL] = assembly_goal_0();

%PRO OTTIMIZZAZIONE
DTOT = VOLFRAC*ones(nelem,1);
dmin   = 1e-6;
dmax   = 1;
DMIN = dmin*ones(nelem,1); 
DMAX = dmax*ones(nelem,1); 

%vincolo di uguaglianza espresso come disuguaglianza
epsi = 10^(-7);
Adi  = 0.25*dx*dy*[ones(1,nelem); -ones(1,nelem)];
Bdi  = [VMAX+epsi; -VMAX+epsi];

%vincolo disuguaglianza
%epsi = 10^(-7);
%Adi  = dx*dy*[ones(1,nelem)];
%Bdi  = [VMAX+eps];

LB = DMIN;
UB = DMAX;

%VINCOLI DI UGUAGLIANZA
Aeq = [];
Beq = [];

%partenza%
load X0
%X0=DTOT/2;
%X0=DMIN;

m = 1+nelem;
n = nelem;
epsimin = 0.0000001;
eeen    = ones(n,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = DTOT;
xold1   = xval;
xold2   = xval;
xmin    = LB;
xmax    = UB;
low     = xmin;
upp     = xmax;
c       = 1000*eeem;
d       = eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1;
kkttol  = 0;

return








%OTTIMIZZA
OPTIONS=optimset('Algorithm','sqp','DerivativeCheck','off','Display','iter','MaxFunEvals',10^8,...
    'GradConstr','off','GradObj','on','TolCon',1e-6,'TolFun',1e-7,...
    'TolX',1e-6,'DiffMinChange',1e-5,'DiffMaxChange', 5e-1,...
   'MaxIter',5,'PlotFcns',@optimplotfval,'UseParallel','always');
%OPTIONS=optimset('Algorithm','sqp','Display','iter','MaxFunEvals',10^8,...
%    'GradConstr','off','GradObj','on','TolCon',1e-6,'TolFun',1e-6,...
%    'TolX',1e-6,'DiffMinChange',1e-6,'DiffMaxChange', 5e-1,...
%    'MaxIter',4,'PlotFcns',@optimplotfval,'UseParallel','always');

%[X,FVAL,EXITFLAG,GRAD,OUTPUT] = fmincon(@ottigoal_1,X0,Adi,Bdi,Aeq,Beq,LB,UB,[],OPTIONS);
[X,FVAL,EXITFLAG,OUTPUT,LAMBDA,GRAD,HESSIAN] = fmincon('ottigoal_2',X0,Adi,Bdi,Aeq,Beq,LB,UB,[],OPTIONS);

X0 = X;
save X0 X0;
%save ALLWORK2

% Deformed mesh
nnod = size(coordinates,1);
Ux = spost(1:2:2*nnod); 
Uy = spost(2:2:2*nnod);
defo = coordinates + [Ux, Uy];


%% PLOT SOLUTION
plotsolution_goal(coordinates,element,X);

plotsolution(coordinates,element,defo,strain,stress);

