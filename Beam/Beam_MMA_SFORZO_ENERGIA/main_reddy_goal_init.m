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
global ndx ndy nnod NUMFREEDISP poisson qexp
global AELE BELE WELE KELE MELE DELE FIXEDDOFS FORCEDDOFS FREEDOFS dx dy VMAX
global A B W M D FORIGINAL iut IOPTION IFILTER ISTRESS pexp HH strain stress spost m n HS SY
global STRESSET NSTRESS ik jk iik jjk iiik jjjk KELE1 KELE2 KELE3 KELE4 KASSEM1

format long

%IOPTION = compliance come: [1 Lext; 2 C_epsi:epsi; 3 D_sigma:epsi
IOPTION = 2;
IFILTER = 1;
ISTRESS = 0; 
NSTRESS = 0; 
pexp = 3;
qexp = 2.9;
SY = 0.035;

%% INPUT
length = 300;                                   % length 
height = 300;                                   % heigth
ndx = 50;                                       % partition in x direction
ndy = 50;                                       % partition in y direction
dx = length/ndx;
dy = height/ndy;
young = 1;                                   % young modulus virgin
youngmin = 1e-8;                              % minimum young modulus
poisson = 0.3;                                % poisson modulus
VOLFRAC = 0.4;
VMAX = length*height*VOLFRAC;

%% GEOMETRY
dx = length/ndx;
dy = height/ndy;
% Mesh 
[coordinates,nnod]=CoordinatesType2(ndx,ndy,dx,dy);
[element,nelem]=ElementType2(ndx,ndy);

if ISTRESS == 1
        %STRESSET = [1:4 4*(ndx-1)+1:4*ndx 4*ndx*(ndy-1)+1:4*ndx*(ndy-1)+4];
        %STRESSET = [1:4 4*ndx*(ndy-1)+1:4*ndx*(ndy-1)+4];
        STRESSET = [1:4];
        NSTRESS = size(STRESSET,2);
end

%
[mc,ngdlu]=CorrispoMC(element,nelem,nnod);  % Corrispondence Matrix displacement 
mc2=CorrispoMC2(element,nelem);             % Corrispondence Matrix strain, stress
iut = [1:4:nelem; 2:4:nelem; 3:4:nelem; 4:4:nelem]';

%dimensionamenti
%% puntatori da 9-9
ik1 = reshape(kron(mc2(1:4:end,:),ones(9,1))',81*nelem/4,1);
ik2 = reshape(kron(mc2(2:4:end,:),ones(9,1))',81*nelem/4,1);
ik3 = reshape(kron(mc2(3:4:end,:),ones(9,1))',81*nelem/4,1);
ik4 = reshape(kron(mc2(4:4:end,:),ones(9,1))',81*nelem/4,1);
jk1 = reshape(kron(mc2(1:4:end,:),ones(1,9))',81*nelem/4,1);
jk2 = reshape(kron(mc2(2:4:end,:),ones(1,9))',81*nelem/4,1);
jk3 = reshape(kron(mc2(3:4:end,:),ones(1,9))',81*nelem/4,1);
jk4 = reshape(kron(mc2(4:4:end,:),ones(1,9))',81*nelem/4,1);
ik = [ik1; ik2; ik3; ik4];
jk = [jk1; jk2; jk3; jk4];
clear ik1 ik2 ik3 ik4 jk1 jk2 jk3 jk4

%% puntatori da 8-9
iik1 = reshape(kron(mc2(1:4:end,:),ones(8,1))',72*nelem/4,1);
iik2 = reshape(kron(mc2(2:4:end,:),ones(8,1))',72*nelem/4,1);
iik3 = reshape(kron(mc2(3:4:end,:),ones(8,1))',72*nelem/4,1);
iik4 = reshape(kron(mc2(4:4:end,:),ones(8,1))',72*nelem/4,1);
jjk1 = reshape(kron(mc(1:4:end,:),ones(1,9))',72*nelem/4,1);
jjk2 = reshape(kron(mc(2:4:end,:),ones(1,9))',72*nelem/4,1);
jjk3 = reshape(kron(mc(3:4:end,:),ones(1,9))',72*nelem/4,1);
jjk4 = reshape(kron(mc(4:4:end,:),ones(1,9))',72*nelem/4,1);
iik = [iik1; iik2; iik3; iik4];
jjk = [jjk1; jjk2; jjk3; jjk4];
clear iik1 iik2 iik3 iik4 jjk1 jjk2 jjk3 jjk4

%% puntatori da 8-8
iiik1 = reshape(kron(mc(1:4:end,:),ones(8,1))',64*nelem/4,1);
iiik2 = reshape(kron(mc(2:4:end,:),ones(8,1))',64*nelem/4,1);
iiik3 = reshape(kron(mc(3:4:end,:),ones(8,1))',64*nelem/4,1);
iiik4 = reshape(kron(mc(4:4:end,:),ones(8,1))',64*nelem/4,1);
jjjk1 = reshape(kron(mc(1:4:end,:),ones(1,8))',64*nelem/4,1);
jjjk2 = reshape(kron(mc(2:4:end,:),ones(1,8))',64*nelem/4,1);
jjjk3 = reshape(kron(mc(3:4:end,:),ones(1,8))',64*nelem/4,1);
jjjk4 = reshape(kron(mc(4:4:end,:),ones(1,8))',64*nelem/4,1);
iiik = [iiik1; iiik2; iiik3; iiik4];
jjjk = [jjjk1; jjjk2; jjjk3; jjjk4];
clear iiik1 iiik2 iiik3 iiik4 jjjk1 jjjk2 jjjk3 jjjk4

ngdls = 3*nnod;
lambda = young*poisson/( (1+poisson)*(1-2*poisson) );
mu = young/(2*(1+poisson));
alpha = 2*mu;
pois(1) = (1-poisson)/((1+poisson)*(1-2*poisson));
pois(2) = poisson/((1+poisson)*(1-2*poisson));
pois(3) = 1/(1+poisson);
pois(4) = 1/(2*(1+poisson));
YOUNG = young*ones(nelem,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GDL LATI
LX1 = 1:2:2*ndx+1; LY1 = LX1+1;
LX2 = 2*ndx+1:2*(ndx+1):2*(ndx+1)*(ndy+1)-1; LY2 = LX2+1;
LX3 = sort(2*(ndx+1)*(ndy+1)-1:-2:2*(ndx+1)*(ndy)+1); LY3 = LX3+1;
LX4 = 1:2*(ndx+1):2*(ndx)*(ndy+1)-3; LY4=LX4+1;

%GDL NODI SPIGOLI
NX1 = min(LX1); NY1 = min(LY1);
NX2 = min(LX2); NY2 = min(LY2); 
NX3 = max(LX3); NY3 = max(LY3);
NX4 = max(LX4); NY4 = max(LY4);

%GDL NODI MEZZERIA
NMX1 = LX1(ndx/2+1); NMY1 = NMX1 + 1;
NMX2 = LX2(ndy/2+1); NMY2 = NMX2 + 1;
NMX3 = LX3(ndx/2+1); NMY3 = NMX3 + 1;
NMX4 = LX4(ndy/2+1); NMY4 = NMX4 + 1;

FIXEDDOFS = [];
%FIXEDDOFS(1,:) = [LX4 NY2]; FIXEDDOFS(2,:) = 0;
FIXEDDOFS(1,:) = [LX4 LY4]; FIXEDDOFS(2,:) = 0;
[nxfixed,nyfixed] = size(FIXEDDOFS);

FORCEDDOFS = [];
FORCEDDOFS(1,:) = NMY2; FORCEDDOFS(2,:) = -1;

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
KELE1 = KELE(:,:,1);
KELE2 = KELE(:,:,2);
KELE3 = KELE(:,:,3);
KELE4 = KELE(:,:,4);

%FILTER
if IFILTER == 1 
    rmin = max(dx,dy);
    [distanze] = distances(coordinates,element);
    HH = zeros(size(distanze));
    HH = max(0,rmin-distanze);
    HS = sum(HH,2);
end

%assemblo una volta per tutte le matrici che non dipendono da Young
[A,B,W,M,D,FORIGINAL] = assembly_goal_00();
KASSEM1 = alpha.*A - alpha.*(B'*D*W + W'*D*B) + W'*D*(alpha.*M)*D*W;
clear iik jjk iiik jjjk
[NUMFREEDISP NUMY] = size(A);
NUMFREEDISP = NUMFREEDISP - nyfixed;

%PRO OTTIMIZZAZIONE
DTOT = VOLFRAC*ones(nelem,1);
dmin   = 1e-9;
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
%load X0
%X0=DTOT/2;
%X0=DMIN;

m = 1+NSTRESS*ISTRESS;
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
c       = 10^4*eeem;
d       = zerom;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 50;
kkttol  = 0.001;
f0valsoglia = 5;

MAINMMA_EXECUTE

return

