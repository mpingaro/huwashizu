% by Marco Pingaro & Paolo Venini

function [A,B,K,M,W,D] = reddy_element(P,lambda,G)

%% Quadrature
%[weight,gauss_x,gauss_y] = gauss_quadrature();
[weight, gauss] = GaussQuad2D(4,4);
gauss_x = gauss(:,1);
gauss_y = gauss(:,2);


%% - Legame
C = [lambda+2*G, 0, 0, lambda; 
    0, 2*G, 0, 0;
    0, 0, 2*G, 0;
    lambda, 0, 0, lambda+2*G];


%% JACOBIAN MATRIX
[DFF,JF] = jacobian(P);

A = zeros(11);    % A
B = zeros(12,11); % B
W = zeros(12,11); % W
K = zeros(12);    % K
M = zeros(12);    % M
D = zeros(12);    % D

%
for i = 1:size(weight,1) % Cycle on gauss points --> Da ottimizzare. (2 con vecchia integrazione
    
   %x = gauss_x(1,i);
   %y = gauss_y(1,i);
   %w = weight(1,i);
   
   x = gauss_x(i,1);
   y = gauss_y(i,1);
   w = weight(i,1);
   
   DFF_i = DFF(:,:,i);
   JF_i  = JF(i);
   
   % Trasformation of the gradient
   % Gradient of shape functions 
   grdu(:,1) = 0.25.*[-(1-y), -(1-x)]*DFF_i; 
   grdu(:,2) = 0.25.*[1-y, -(1+x)]*DFF_i;
   grdu(:,3) = 0.25.*[1+y, 1+x]*DFF_i;
   grdu(:,4) = 0.25.*[-(1+y), 1-x]*DFF_i;
   % Grad of Boubble functions
   % First Boubble function
   grdu(:,5) = [-2*x*(1-y^2), -2*y*(1-x^2)]*DFF_i;
   % Mixed Boubble function
   grdu(:,6) = [-2*x*(y-y^3-1+y^2), (1-3*y^2+2*y)*(1-x^2)]*DFF_i;
   grdu(:,7) = [(1-3*x^2+2*x)*(1-y^2), -2*y*(x-x^3-1+x^2)]*DFF_i;
   
   
   epsi = [grdu(1,1), 0, grdu(1,2), 0, grdu(1,3), 0, grdu(1,4), 0, grdu(1,5), 0, grdu(1,6);
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2,...
           grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2, grdu(2,5)/2, grdu(1,5)/2, (grdu(1,7)+grdu(2,6))/2;
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2,...
           grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2, grdu(2,5)/2, grdu(1,5)/2, (grdu(1,7)+grdu(2,6))/2;
           0, grdu(2,1), 0, grdu(2,2), 0, grdu(2,3), 0, grdu(2,4), 0, grdu(2,5), grdu(2,7)];
   
       
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
