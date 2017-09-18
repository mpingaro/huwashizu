% by Marco Pingaro & Paolo Venini

function load = body_load_error(P,l,m)

%% Quadrature
%[weight,gauss_x,gauss_y] = gauss_quadrature();
[weight, gauss] = GaussQuad2D(4,4);
gauss_x = gauss(:,1);
gauss_y = gauss(:,2);


%% JACOBIAN MATRIX
[DFF,JF] = jacobian(P);

load = zeros(10,1);
for i = 1:size(weight,1) % 2 prima
    
   [x,y] = map_quad(P,gauss(i,:));
   xi = gauss_x(i,1);
   eta = gauss_y(i,1);
   w = weight(i,1);
   
   % Shape Functions
   psi(1) = 0.25*(1-xi)*(1-eta);
   psi(2) = 0.25*(1+xi)*(1-eta);
   psi(3) = 0.25*(1+xi)*(1+eta);
   psi(4) = 0.25*(1-xi)*(1+eta);
   %psi(5) = 4*psi(1)*(1-xi^2)*(1-eta^2);   
   psi(5) = (1+xi+eta)*(1-xi^2)*(1-eta^2);

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
   
   A = 2/(1+l);
   %B = 0.5*A*sin(pi*x)*sin(pi*y);
   b = 1/25;

   g(1,1) = b*(pi^2*(4*sin(2*pi*y)*(-1+2*cos(2*pi*x))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 
   g(1,2) = b*(pi^2*(4*sin(2*pi*x)*( 1-2*cos(2*pi*y))-cos(pi*(x+y))+A*sin(pi*x)*sin(pi*y))); 

   for j = 1:10
       load(j,1) = load(j,1) + w*( g(1,1)*v(j,1)+g(1,2)*v(j,2) )*JF(i); 
   end
   
end
