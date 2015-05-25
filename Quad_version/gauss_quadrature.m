% by Marco Pingaro & Paolo Venini

function [weight,gauss_x,gauss_y] = gauss_quadrature()

%% QUADRATURA DI GAUSS (grado di precisione 5)
% Weight of quadrature
weight(1,1)=0.112500000000000; weight(1,2)=0.066197076394253; 
weight(1,3)=0.066197076394253; weight(1,4)=0.066197076394253; 
weight(1,5)=0.062969590272414; weight(1,6)=0.062969590272414; 
weight(1,7)=0.062969590272414;
% Point value 
gauss_x(1,1) = 0.333333333333333; gauss_x(1,2) = 0.470142064105115; 
gauss_x(1,3) = 0.059715871789770; gauss_x(1,4) = 0.470142064105115; 
gauss_x(1,5) = 0.101286507323456; gauss_x(1,6) = 0.797426985353087;
gauss_x(1,7) = 0.101286507323456; 
gauss_y(1,1) = 0.333333333333333; gauss_y(1,2) = 0.470142064105115; 
gauss_y(1,3) = 0.470142064105115; gauss_y(1,4) = 0.059715871789770;
gauss_y(1,5) = 0.101286507323456; gauss_y(1,6) = 0.101286507323456;
gauss_y(1,7) = 0.797426985353087;

end
