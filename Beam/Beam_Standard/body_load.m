% by Marco Pingaro & Paolo Venini

function load = body_load(P,g)

%% Quadrature
[weight,gauss_x,gauss_y] = gauss_quadrature();

%% JACOBIAN MATRIX
DF(1,1) = P(3)-P(1);
DF(1,2) = P(5)-P(1);
DF(2,1) = P(4)-P(2);
DF(2,2) = P(6)-P(2);
% Determinant of Jacobian matrix
JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);

load = zeros(8,1);
for i = 1:size(weight,2)
    
   x = gauss_x(1,i);
   y = gauss_y(1,i);
   w = weight(1,i);
   
   psi(1) = 1-x-y;
   psi(2) = x;
   psi(3) = y;
   psi(4) = 27*psi(1)*psi(2)*psi(3);
   
   v(1,:) = [psi(1), 0];
   v(2,:) = [0, psi(1)];
   v(3,:) = [psi(2), 0];
   v(4,:) = [0, psi(2)];
   v(5,:) = [psi(3), 0];
   v(6,:) = [0, psi(3)];
   v(7,:) = [psi(4), 0];
   v(8,:) = [0, psi(4)];   
   
   for j = 1:8
       load(j,1) = load(j,1) + w*( g(1,1)*v(j,1)+g(1,2)*v(j,2) )*JF; 
   end
   
end