% by Marco Pingaro & Paolo Venini

function load = body_load(P,g)

%% Quadrature
[weight,gauss_x,gauss_y] = gauss_quadrature();

%% JACOBIAN MATRIX
[DFF,JF] = jacobian(P);

load = zeros(12,1);
for i = 1:size(weight,2)
    
   x = gauss_x(1,i);
   y = gauss_y(1,i);
   w = weight(1,i);
   
   % Shape Functions
   psi(1) = 0.25*(1-x)*(1-y);
   psi(2) = 0.25*(1+x)*(1-y);
   psi(3) = 0.25*(1+x)*(1+y);
   psi(4) = 0.25*(1-x)*(1+y);
   psi(5) = (1-x^2)*(1-y^2);
   psi(6) = (x+y)*psi(5) ;
    
   v(1,:)  = [psi(1), 0];
   v(2,:)  = [0, psi(1)];
   v(3,:)  = [psi(2), 0];
   v(4,:)  = [0, psi(2)];
   v(5,:)  = [psi(3), 0];
   v(6,:)  = [0, psi(3)];
   v(7,:)  = [psi(4), 0];
   v(8,:)  = [0, psi(4)];
   v(9,:)  = [psi(5), 0];
   v(10,:) = [0, psi(5)];
   v(11,:) = [psi(6),0];
   v(12,:) = [0,psi(6)];
   
   for j = 1:12
       load(j,1) = load(j,1) + w*( g(1,1)*v(j,1)+g(1,2)*v(j,2) )*JF(i); 
   end
   
end