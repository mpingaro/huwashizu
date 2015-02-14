% by Marco Pingaro & Paolo Venini

function [A,B,K,M,W,D] = reddy_3Delement(point,lambda,G)

%% Quadrature
[xa,ya,za,wt] = TetQuadDat(3);
%[xa,ya,za,wt]=tetraquad(5,[0,0,0;1,0,0;0,1,0;0,0,1]);

%% - Legame
Cblk1 = [lambda+2*G, lambda, lambda;
      lambda, lambda+2*G, lambda;
      lambda, lambda, lambda+2*G];
Cblk2 = diag([2*G,2*G,2*G,2*G,2*G,2*G],0);
C = blkdiag(Cblk1,Cblk2);

%% Jacobian matrix 
[jac,inv_jac,djac] = jacobian_tetra(point);

A = zeros(15,15); % A
B = zeros(24,15); % B
W = zeros(24,15); % W
K = zeros(24,24); % K
M = zeros(24,24); % M
D = zeros(24,24); % D

% Trasformation of the gradient functions (constant part)
grd(:,1) = inv_jac*[-1; -1; -1];
grd(:,2) = inv_jac*[ 1;  0;  0];
grd(:,3) = inv_jac*[ 0;  1;  0];
grd(:,4) = inv_jac*[ 0;  0;  1];

for i = 1:size(wt,2) % Cycle on gauss points --> Da ottimizzare.
    
   x = xa(i); y = ya(i); z = za(i); w = wt(i);
   
   % Grad of Boubble function
   grd(:,5) = inv_jac*[256*y*z*(1-2*x-y-z); 256*x*z*(1-x-2*y-z); 256*x*y*(1-x-y-2*z)];   
   
   % first node
   epsi(:,1) = [grd(1,1), 0, 0, grd(2,1)/2, grd(2,1)/2, grd(3,1)/2, grd(3,1)/2, 0, 0];
   epsi(:,2) = [0, grd(2,1), 0, grd(1,1)/2, grd(1,1)/2, 0, 0, grd(3,1)/2, grd(3,1)/2];
   epsi(:,3) = [0, 0, grd(3,1), 0, 0, grd(1,1)/2, grd(1,1)/2, grd(2,1)/2, grd(2,1)/2];
   % second node
   epsi(:,4) = [grd(1,2), 0, 0, grd(2,2)/2, grd(2,2)/2, grd(3,2)/2, grd(3,2)/2, 0, 0];
   epsi(:,5) = [0, grd(2,2), 0, grd(1,2)/2, grd(1,2)/2, 0, 0, grd(3,2)/2, grd(3,2)/2];
   epsi(:,6) = [0, 0, grd(3,2), 0, 0, grd(1,2)/2, grd(1,2)/2, grd(2,2)/2, grd(2,2)/2];
   % third node
   epsi(:,7) = [grd(1,3), 0, 0, grd(2,3)/2, grd(2,3)/2, grd(3,3)/2, grd(3,3)/2, 0, 0];
   epsi(:,8) = [0, grd(2,3), 0, grd(1,3)/2, grd(1,3)/2, 0, 0, grd(3,3)/2, grd(3,3)/2];
   epsi(:,9) = [0, 0, grd(3,3), 0, 0, grd(1,3)/2, grd(1,3)/2, grd(2,3)/2, grd(2,3)/2];
   % fourth node
   epsi(:,10) = [grd(1,4), 0, 0, grd(2,4)/2, grd(2,4)/2, grd(3,4)/2, grd(3,4)/2, 0, 0];
   epsi(:,11) = [0, grd(2,4), 0, grd(1,4)/2, grd(1,4)/2, 0, 0, grd(3,4)/2, grd(3,4)/2];
   epsi(:,12) = [0, 0, grd(3,4), 0, 0, grd(1,4)/2, grd(1,4)/2, grd(2,4)/2, grd(2,4)/2];
   % Bouble functions
   epsi(:,13) = [grd(1,5), 0, 0, grd(2,5)/2, grd(2,5)/2, grd(3,5)/2, grd(3,5)/2, 0, 0];
   epsi(:,14) = [0, grd(2,5), 0, grd(1,5)/2, grd(1,5)/2, 0, 0, grd(3,5)/2, grd(3,5)/2];
   epsi(:,15) = [0, 0, grd(3,5), 0, 0, grd(1,5)/2, grd(1,5)/2, grd(2,5)/2, grd(2,5)/2];
    
   %% Shape function strain
   d(1) = 1-x-y-z;
   d(2) = x;
   d(3) = y;
   d(4) = z;
   
   % d_xx, d_yy, d_zz, d_xy, d_yx, d_xz, d_zx, d_yz, d_zy
   % first node
   strain(:,1) = [d(1), 0, 0, 0, 0, 0, 0, 0, 0];
   strain(:,2) = [0, d(1), 0, 0, 0, 0, 0, 0, 0];
   strain(:,3) = [0, 0, d(1), 0, 0, 0, 0, 0, 0];
   strain(:,4) = [0, 0, 0, d(1), d(1), 0, 0, 0, 0];
   strain(:,5) = [0, 0, 0, 0, 0, d(1), d(1), 0, 0];
   strain(:,6) = [0, 0, 0, 0, 0, 0, 0, d(1), d(1)];
   % second node
   strain(:,7) = [d(2), 0, 0, 0, 0, 0, 0, 0, 0];
   strain(:,8) = [0, d(2), 0, 0, 0, 0, 0, 0, 0];
   strain(:,9) = [0, 0, d(2), 0, 0, 0, 0, 0, 0];
   strain(:,10) = [0, 0, 0, d(2), d(2), 0, 0, 0, 0];
   strain(:,11) = [0, 0, 0, 0, 0, d(2), d(2), 0, 0];
   strain(:,12) = [0, 0, 0, 0, 0, 0, 0, d(2), d(2)];
   % third node
   strain(:,13) = [d(3), 0, 0, 0, 0, 0, 0, 0, 0];
   strain(:,14) = [0, d(3), 0, 0, 0, 0, 0, 0, 0];
   strain(:,15) = [0, 0, d(3), 0, 0, 0, 0, 0, 0];
   strain(:,16) = [0, 0, 0, d(3), d(3), 0, 0, 0, 0];
   strain(:,17) = [0, 0, 0, 0, 0, d(3), d(3), 0, 0];
   strain(:,18) = [0, 0, 0, 0, 0, 0, 0, d(3), d(3)];
   % fourth node
   strain(:,19) = [d(4), 0, 0, 0, 0, 0, 0, 0, 0];
   strain(:,20) = [0, d(4), 0, 0, 0, 0, 0, 0, 0];
   strain(:,21) = [0, 0, d(4), 0, 0, 0, 0, 0, 0];
   strain(:,22) = [0, 0, 0, d(4), d(4), 0, 0, 0, 0];
   strain(:,23) = [0, 0, 0, 0, 0, d(4), d(4), 0, 0];
   strain(:,24) = [0, 0, 0, 0, 0, 0, 0, d(4), d(4)];
   
   %% Shape function stress
   t(1) = 4-5*x-5*y-5*z;
   t(2) = 5*x-1;
   t(3) = 5*y-1;
   t(4) = 5*z-1;
   
   % t_xx, t_yy, t_zz, t_xy, t_yx, t_xz, t_zx, t_yz, t_zy
   % first node
   tau(:,1) = [t(1), 0, 0, 0, 0, 0, 0, 0, 0];
   tau(:,2) = [0, t(1), 0, 0, 0, 0, 0, 0, 0];
   tau(:,3) = [0, 0, t(1), 0, 0, 0, 0, 0, 0];
   tau(:,4) = [0, 0, 0, t(1), t(1), 0, 0, 0, 0];
   tau(:,5) = [0, 0, 0, 0, 0, t(1), t(1), 0, 0];
   tau(:,6) = [0, 0, 0, 0, 0, 0, 0, t(1), t(1)];
   % second node
   tau(:,7) = [t(2), 0, 0, 0, 0, 0, 0, 0, 0];
   tau(:,8) = [0, t(2), 0, 0, 0, 0, 0, 0, 0];
   tau(:,9) = [0, 0, t(2), 0, 0, 0, 0, 0, 0];
   tau(:,10) = [0, 0, 0, t(2), t(2), 0, 0, 0, 0];
   tau(:,11) = [0, 0, 0, 0, 0, t(2), t(2), 0, 0];
   tau(:,12) = [0, 0, 0, 0, 0, 0, 0, t(2), t(2)];
   % third node
   tau(:,13) = [t(3), 0, 0, 0, 0, 0, 0, 0, 0];
   tau(:,14) = [0, t(3), 0, 0, 0, 0, 0, 0, 0];
   tau(:,15) = [0, 0, t(3), 0, 0, 0, 0, 0, 0];
   tau(:,16) = [0, 0, 0, t(3), t(3), 0, 0, 0, 0];
   tau(:,17) = [0, 0, 0, 0, 0, t(3), t(3), 0, 0];
   tau(:,18) = [0, 0, 0, 0, 0, 0, 0, t(3), t(3)];
   % fourth node
   tau(:,19) = [t(4), 0, 0, 0, 0, 0, 0, 0, 0];
   tau(:,20) = [0, t(4), 0, 0, 0, 0, 0, 0, 0];
   tau(:,21) = [0, 0, t(4), 0, 0, 0, 0, 0, 0];
   tau(:,22) = [0, 0, 0, t(4), t(4), 0, 0, 0, 0];
   tau(:,23) = [0, 0, 0, 0, 0, t(4), t(4), 0, 0];
   tau(:,24) = [0, 0, 0, 0, 0, 0, 0, t(4), t(4)];
   
   %% Matrix:  A( epsi(u),apsi(v) ) Stored full matrix
   A = A + w.*epsi'*epsi.*djac;
   %% Matrix:  B( e(u),apsi(v) ) Stored full matrix
   B = B + w.*strain'*epsi.*djac;
   %% Matrix:  W( tau(u),apsi(v) ) Stored full matrix
   W = W + w.*tau'*epsi.*djac;
   %% Matrix:  K( C*e(u),d(v) ) Stored full matrix (Symmetric)
   K = K + w.*(C*strain)'*strain.*djac;
   %% Matrix:  M( e(u),d(v) ) Stored full matrix (Symmetric)
   M = M + w.*strain'*strain.*djac;
   %% Matrix:  D( tau(u),d(v) ) Stored full matrix
   D = D + w.*tau'*strain.*djac;
end

end