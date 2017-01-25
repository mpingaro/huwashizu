% by Marco Pingaro & Paolo Venini

function M = mass_matrix(P,rho)

%% Quadrature
[weight,gauss_x,gauss_y] = gauss_quadrature();

%% JACOBIAN MATRIX
%% JACOBIAN MATRIX
[DFF,JF] = jacobian(P);

M = zeros(10);          % Elementar mass matrix
for i = 1:size(weight,2)   % Cycle on gauss points --> Da ottimizzare.
   x = gauss_x(1,i);
   y = gauss_y(1,i);
   w = weight(1,i);
   % Shape Functions
   psi(1) = 0.25*(1-x)*(1-y);
   psi(2) = 0.25*(1+x)*(1-y);
   psi(3) = 0.25*(1+x)*(1+y);
   psi(4) = 0.25*(1-x)*(1+y);
   psi(5) = (1-x^2)*(1-y^2);
   %
   u = [psi(1), 0, psi(2), 0, psi(3), 0, psi(4), 0, psi(5), 0; 
       0, psi(1), 0, psi(2), 0, psi(3), 0, psi(4), 0, psi(5)];
   
   %% Matrix:  M( rho* u,v ) Stored full matrix
   M = M + w.*( rho.*u'*u ).*JF(i);
end

end
