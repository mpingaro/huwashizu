% by Marco Pingaro & Paolo Venini

function load = body_load(point,g)

%% Quadrature
[xa,ya,za,wt] = TetQuadDat(4); % Max order of quadrature

%% Jacobian matrix 
[jac,inv_jac,djac] = jacobian_tetra(point);

load = zeros(15,1);
for i = 1:size(wt,2)
    
   x = xa(i); y = ya(i); z = za(i); w = wt(i);
   
   %% Shape function Displacement
   psi(1) = 1-x-y-z;
   psi(2) = x;
   psi(3) = y;
   psi(4) = z;
   psi(5) = 256*x*y*z*(1-x-y-z);
   
   v(1,:) = [psi(1), 0, 0];
   v(2,:) = [0, psi(1), 0];
   v(3,:) = [0, 0, psi(1)];
   v(4,:) = [psi(2), 0, 0];
   v(5,:) = [0, psi(2), 0];
   v(6,:) = [0, 0, psi(2)];
   v(7,:) = [psi(3), 0, 0];
   v(8,:) = [0, psi(3), 0];
   v(9,:) = [0, 0, psi(3)];
   v(10,:) = [psi(4), 0, 0];
   v(11,:) = [0, psi(4), 0];
   v(12,:) = [0, 0, psi(4)];
   v(13,:) = [psi(5), 0, 0];
   v(14,:) = [0, psi(5), 0];
   v(15,:) = [0, 0, psi(5)];
   
   for j = 1:15
       load(j,1) = load(j,1) + w*( g(1,1)*v(j,1)+g(1,2)*v(j,2)+g(1,3)*v(j,3) )*djac; 
   end
   
end