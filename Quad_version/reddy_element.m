% by Marco Pingaro & Paolo Venini

function [A,B,K,M,W,D] = reddy_element(P,lambda,G)

%% Quadrature
[weight,gauss_x,gauss_y] = gauss_quadrature();

%% - Legame
C = [lambda+2*G, 0, 0, lambda; 
    0, 2*G, 0, 0;
    0, 0, 2*G, 0;
    lambda, 0, 0, lambda+2*G];

%% JACOBIAN MATRIX
DF(1,1) = P(3)-P(1);
DF(1,2) = P(5)-P(1);
DF(2,1) = P(4)-P(2);
DF(2,2) = P(6)-P(2);
% Determinant of Jacobian matrix
JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);
% Inverse transpose of Jacobian matrix
DFF(1,1) = DF(2,2)/JF;
DFF(1,2) = -DF(2,1)/JF;
DFF(2,1) = -DF(1,2)/JF;
DFF(2,2) = DF(1,1)/JF;

A = zeros(8,8); % A
B = zeros(9,8); % B
W = zeros(9,8); % W
K = zeros(9,9); % K
M = zeros(9,9); % M
D = zeros(9,9); % D

% Trasformation of the gradient (constant part)
grdu(:,1) = DFF*[-1;-1]; 
grdu(:,2) = DFF*[1; 0];
grdu(:,3) = DFF*[0; 1];
%
for i = 1:size(weight,2) % Cycle on gauss points --> Da ottimizzare.
    
   x = gauss_x(1,i);
   y = gauss_y(1,i);
   w = weight(1,i);
   
   % Grad of Boubble function
   grdu(:,4) = DFF*[y-2*x*y-y*y; x-2*x*y-x*x]*27;
   
   epsi = [grdu(1,1), 0, grdu(1,2), 0, grdu(1,3), 0, grdu(1,4), 0;
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2, grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2;
           grdu(2,1)/2, grdu(1,1)/2, grdu(2,2)/2, grdu(1,2)/2, grdu(2,3)/2, grdu(1,3)/2, grdu(2,4)/2, grdu(1,4)/2;
           0, grdu(2,1), 0, grdu(2,2), 0, grdu(2,3), 0, grdu(2,4)];
       
   d(1,1) = 1-x-y;
   d(1,2) = x;
   d(1,3) = y;
   
   strain = [d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0, 0;
             0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0; 
             0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3), 0;
             0, 0, d(1,1), 0, 0, d(1,2), 0, 0, d(1,3)];
         
   t(1,1) = 3-4*x-4*y; 
   t(1,2) = 4*x-1; 
   t(1,3) = 4*y-1;
   
   tau = [t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0, 0;
          0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0; 
          0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3), 0;
          0, 0, t(1,1), 0, 0, t(1,2), 0, 0, t(1,3)];
   
   %% Matrix:  A( epsi(u),apsi(v) ) Stored full matrix
   A = A + w.*epsi'*epsi.*JF;
   %% Matrix:  B( e(u),apsi(v) ) Stored full matrix
   B = B + w.*strain'*epsi.*JF;
   %% Matrix:  W( tau(u),apsi(v) ) Stored full matrix
   W = W + w.*tau'*epsi.*JF;
   %% Matrix:  K( C*e(u),d(v) ) Stored full matrix (Symmetric)
   K = K + w.*(C*strain)'*strain.*JF;
   %% Matrix:  M( e(u),d(v) ) Stored full matrix (Symmetric)
   M = M + w.*strain'*strain.*JF;
   %% Matrix:  D( tau(u),d(v) ) Stored full matrix
   D = D + w.*tau'*strain.*JF;
end

end



% OLD VERSION
% %% Quadrature
% [weight,gauss_x,gauss_y] = gauss_quadrature();
% 
% %% - Legame
% C = [lambda+2*G, lambda, 2*G];
% 
% %% JACOBIAN MATRIX
% DF(1,1) = P(3)-P(1);
% DF(1,2) = P(5)-P(1);
% DF(2,1) = P(4)-P(2);
% DF(2,2) = P(6)-P(2);
% % Determinant of Jacobian matrix
% JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);
% % Inverse transpose of Jacobian matrix
% DFF(1,1) = DF(2,2)/JF;
% DFF(1,2) = -DF(2,1)/JF;
% DFF(2,1) = -DF(1,2)/JF;
% DFF(2,2) = DF(1,1)/JF;
% 
% A = zeros(8,8); % A
% B = zeros(9,8); % B
% W = zeros(9,8); % W
% K = zeros(9,9); % K
% M = zeros(9,9); % M
% D = zeros(9,9); % D
% 
% % Trasformation of the gradient (constant part)
% grdu(:,1) = DFF*[-1;-1]; 
% grdu(:,2) = DFF*[1; 0];
% grdu(:,3) = DFF*[0; 1];
% %
% for i = 1:size(weight,2) % Cycle on gauss points --> Da ottimizzare.
%     
%    x = gauss_x(1,i);
%    y = gauss_y(1,i);
%    w = weight(1,i);
%    
%    % Grad of Boubble function
%    grdu(:,4) = DFF*[y-2*x*y-y*y; x-2*x*y-x*x]*27;
%
%    %% Matrix:  A( epsi(u),apsi(v) ) Stored full matrix
%    % Shape functions epsilon(u)
%    epsi(1,:) = [grdu(1,1), 0, grdu(2,1)/2];
%    epsi(2,:) = [0, grdu(2,1), grdu(1,1)/2];
%    epsi(3,:) = [grdu(1,2), 0, grdu(2,2)/2];
%    epsi(4,:) = [0, grdu(2,2), grdu(1,2)/2];
%    epsi(5,:) = [grdu(1,3), 0, grdu(2,3)/2];
%    epsi(6,:) = [0, grdu(2,3), grdu(1,3)/2];
%    epsi(7,:) = [grdu(1,4), 0, grdu(2,4)/2];
%    epsi(8,:) = [0, grdu(2,4), grdu(1,4)/2];   
%
%    for nb_i=1:8
%        for nb_j=1:8
%                A(nb_i,nb_j) = A(nb_i,nb_j) + w*( epsi(nb_i,1)*epsi(nb_j,1)+...
%                    epsi(nb_i,2)*epsi(nb_j,2) + 2*epsi(nb_i,3)*epsi(nb_j,3))*JF;
%        end
%    end
% 
%    %% Matrix:  B( e(u),apsi(v) ) Stored full matrix
%    % Shape functions strains
%    d(1,1) = 1-x-y;
%    d(1,2) = x;
%    d(1,3) = y;
%    strain(1,:) = [d(1,1),0,0];
%    strain(2,:) = [0,0,d(1,1)];
%    strain(3,:) = [0,d(1,1),0];
%    strain(4,:) = [d(1,2),0,0];
%    strain(5,:) = [0,0,d(1,2)];
%    strain(6,:) = [0,d(1,2),0];
%    strain(7,:) = [d(1,3),0,0];
%    strain(8,:) = [0,0,d(1,3)];
%    strain(9,:) = [0,d(1,3),0];
%    for nb_i=1:9
%        for nb_j=1:8
%                B(nb_i,nb_j) = B(nb_i,nb_j) + w*( strain(nb_i,1)*epsi(nb_j,1)+...
%                    strain(nb_i,2)*epsi(nb_j,2)+2*strain(nb_i,3)*epsi(nb_j,3))*JF;
%        end
%    end
%   
%    %% Matrix:  W( tau(u),apsi(v) ) Stored full matrix
%    % Shape functions stress
%    t(1,1) = 3-4*x-4*y; 
%    t(1,2) = 4*x-1; 
%    t(1,3) = 4*y-1;
%    tau(1,:) = [t(1,1),0,0];
%    tau(2,:) = [0,0,t(1,1)];
%    tau(3,:) = [0,t(1,1),0];
%    tau(4,:) = [t(1,2),0,0];
%    tau(5,:) = [0,0,t(1,2)];
%    tau(6,:) = [0,t(1,2),0];
%    tau(7,:) = [t(1,3),0,0];
%    tau(8,:) = [0,0,t(1,3)];
%    tau(9,:) = [0,t(1,3),0];
%    
%    for nb_i=1:9
%        for nb_j=1:8
%                W(nb_i,nb_j) = W(nb_i,nb_j) + w*( tau(nb_i,1)*epsi(nb_j,1)+...
%                    tau(nb_i,2)*epsi(nb_j,2)+2*tau(nb_i,3)*epsi(nb_j,3))*JF;
%        end
%    end
%    
%    for nb_i=1:9
%        for nb_j=1:9
%            %% Matrix:  K( C*e(u),d(v) ) Stored full matrix (Symmetric)
%            K(nb_i,nb_j) = K(nb_i,nb_j) +...
%                w*( (C(1,1)*strain(nb_i,1)+C(1,2)*strain(nb_i,2))*strain(nb_j,1)+...
%                (C(1,1)*strain(nb_i,2)+C(1,2)*strain(nb_i,1))*strain(nb_j,2)+...
%                2*( C(1,3)*strain(nb_i,3) )*strain(nb_j,3) )*JF;
%            
%            %% Matrix:  M( e(u),d(v) ) Stored full matrix (Symmetric)
%            M(nb_i,nb_j) = M(nb_i,nb_j) + w*( strain(nb_i,1)*strain(nb_j,1)+...
%                strain(nb_i,2)*strain(nb_j,2)+2*strain(nb_i,3)*strain(nb_j,3))*JF;
%        end
%        %% Matrix:  D( tau(u),d(v) ) Stored full matrix
%        D(nb_i,nb_i) = D(nb_i,nb_i) + w*( tau(nb_i,1)*strain(nb_i,1)+...
%            tau(nb_i,2)*strain(nb_i,2)+2*tau(nb_i,3)*strain(nb_i,3))*JF;
%    end
% end