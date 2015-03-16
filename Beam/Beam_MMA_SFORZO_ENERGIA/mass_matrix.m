% by Marco Pingaro & Paolo Venini

function M = mass_matrix(P,rho)

%% Quadrature
[weight,gauss_x,gauss_y] = gauss_quadrature();

%% JACOBIAN MATRIX
DF(1,1) = P(3)-P(1);
DF(1,2) = P(5)-P(1);
DF(2,1) = P(4)-P(2);
DF(2,2) = P(6)-P(2);
% Determinant of Jacobian matrix
JF = DF(1,1)*DF(2,2)-DF(1,2)*DF(2,1);

M = zeros(8,8);          % Elementar mass matrix
for i = 1:size(weight,2) % Cycle on gauss points --> Da ottimizzare.
   x = gauss_x(1,i);
   y = gauss_y(1,i);
   w = weight(1,i);
   % Grad of Boubble function
   psi(1) = 1-x-y;
   psi(2) = x;
   psi(3) = y;
   psi(4) = 27*psi(1)*psi(1)*psi(1);
   %
   u = [psi(1), 0, psi(2), 0, psi(3), 0, psi(4), 0; 
       0, psi(1), 0, psi(2), 0, psi(3), 0, psi(4)];
   
   %% Matrix:  M( rho* u,v ) Stored full matrix
   M = M + w.*( rho.*u'*u ).*JF;
end

end