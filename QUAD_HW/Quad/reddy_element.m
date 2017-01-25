% by Marco Pingaro & Paolo Venini

function [A,B,K,M,W,D] = reddy_element(P,lambda,G)

%% Quadrature
%[weight,gauss_x,gauss_y] = gauss_quadrature();
[weight, gauss] = GaussQuad2D(9,9);
gauss_x = gauss(:,1);
gauss_y = gauss(:,2);


%% - Legame
C = [lambda+2*G, 0, 0, lambda; 
    0, 2*G, 0, 0;
    0, 0, 2*G, 0;
    lambda, 0, 0, lambda+2*G];


%% JACOBIAN MATRIX
[DFF,JF] = jacobian(P);

A = zeros(8);    % A
B = zeros(12,8); % B
W = zeros(12,8); % W
K = zeros(12);    % K
M = zeros(12);    % M
D = zeros(12);    % D
%
for i = 1:size(weight,1) % Cycle on gauss points --> Da ottimizzare. (2 prima)
 
   x = gauss_x(i,1);
   y = gauss_y(i,1);
   w = weight(i,1);
   
   DFF_i = DFF(:,:,i);
   JF_i  = JF(i);
   
   % Trasformation of the gradient
   % Gradient of shape functions 
   grdu(:,1) = DFF_i*[-(1-y); -(1-x)].*0.25; 
   grdu(:,2) = DFF_i*[1-y; -(1+x)].*0.25;
   grdu(:,3) = DFF_i*[1+y; 1+x].*0.25;
   grdu(:,4) = DFF_i*[-(1+y); 1-x].*0.25;

   epsi = [grdu(1,1), 0, grdu(1,2), 0, grdu(1,3), 0, grdu(1,4), 0;
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2,...
           grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2;
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2,...
           grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2;
           0, grdu(2,1), 0, grdu(2,2), 0, grdu(2,3), 0, grdu(2,4)];
    
   d(1,1) = 0.25*(1-x)*(1-y) ;
   d(1,2) = 0.25*(1+x)*(1-y) ;
   d(1,3) = 0.25*(1+x)*(1+y) ;
   d(1,4) = 0.25*(1-x)*(1+y) ;
   
   strain = [d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0, 0, d(1,4), 0, 0;
             0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0, 0, d(1,4), 0; 
             0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0, 0, d(1,4), 0;
             0, 0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0, 0, d(1,4)];
         
   t(1,1) = 1 - 3*x - 3*y + 9*x*y;
   t(1,2) = 1 + 3*x - 3*y - 9*x*y;
   t(1,3) = 1 + 3*x + 3*y + 9*x*y;
   t(1,4) = 1 - 3*x + 3*y + 9*x*y;
   
   tau = [t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0, 0, t(1,4), 0, 0;
          0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0, 0, t(1,4), 0; 
          0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0, 0, t(1,4), 0;
          0, 0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0, 0, t(1,4)];
   
   %% Matrix:  A( epsi(u),apsi(v) ) Stored full matrix
   A = A + w.*epsi'*epsi.*JF_i;
   %% Matrix:  B( e(u),apsi(v) ) Stored full matrix
   B = B + w.*strain'*epsi.*JF_i;
   %% Matrix:  W( tau(u),apsi(v) ) Stored full matrix
   W = W + w.*tau'*epsi.*JF_i;
   %% Matrix:  K( C*e(u),d(v) ) Stored full matrix (Symmetric)
   K = K + w.*(C*strain)'*strain.*JF_i;
   %% Matrix:  M( e(u),d(v) ) Stored full matrix (Symmetric)
   M = M + w.*strain'*strain.*JF_i;
   %% Matrix:  D( tau(u),d(v) ) Stored full matrix
   D = D + w.*tau'*strain.*JF_i;
end

end
