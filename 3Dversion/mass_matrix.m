% by Marco Pingaro & Paolo Venini

function M = mass_matrix(P,rho)

%% QUADRATURE
[xa,ya,za,wt] = TetQuadDat(4);

%% JACOBIAN MATRIX
[jac,inv_jac,djac] = jacobian_tetra(P);

M = zeros(15,15);          % Elementar mass matrix
for i = 1:size(wt,2)       % Cycle on gauss points --> Da ottimizzare.
   x = xa(i); y = ya(i); z = za(i); w = wt(i);
   % SHAPE FUNCTION
   psi(1) = 1-x-y-z;
   psi(2) = x;
   psi(3) = y;
   psi(4) = z;
   psi(5) = 256*x*y*z*(1-x-y-z);
   %
   u = [psi(1), 0, 0, psi(2), 0, 0, psi(3), 0, 0, psi(4), 0, 0, psi(5), 0, 0; 
        0, psi(1), 0, 0, psi(2), 0, 0, psi(3), 0, 0, psi(4), 0, 0, psi(5), 0;
        0, 0, psi(1), 0, 0, psi(2), 0, 0, psi(3), 0, 0, psi(4), 0, 0, psi(5)];
   
   %% Matrix:  M( rho* u,v ) Stored full matrix
   M = M + w.*( rho.*u'*u ).*djac;
end

end