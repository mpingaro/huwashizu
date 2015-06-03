% by Marco Pingaro & Paolo Venini

function [weight,gauss_x,gauss_y] = gauss_quadrature_9()

 
%% QUADRATURA DI GAUSS
% Point value
% Direction x
gauss_x(1,1) = -0.774596669241483 ; gauss_x(1,2) = -0.774596669241483 ;
gauss_x(1,3) = -0.774596669241483 ; gauss_x(1,4) = 0.0 ; 
gauss_x(1,5) = 0.0 ; gauss_x(1,6) = 0.0 ; gauss_x(1,7) = 0.774596669241483 ; 
gauss_x(1,8) = 0.774596669241483 ; gauss_x(1,9) = 0.774596669241483 ;
% Direction y
gauss_y(1,1) = -0.774596669241483 ; gauss_y(1,2) = 0.0 ; 
gauss_y(1,3) = 0.774596669241483 ; gauss_y(1,4) = -0.774596669241483 ; 
gauss_y(1,5) = 0 ; gauss_y(1,6) = 0.774596669241483 ; 
gauss_y(1,7) = -0.774596669241483 ; gauss_y(1,8) = 0.0 ; 
gauss_y(1,9) = 0.774596669241483 ;
% Weight of Quadrature
weight(1,1) = 0.308641975308642 ; weight(1,2) = 0.493827160493827 ;
weight(1,3) = 0.308641975308642 ; weight(1,4) = 0.493827160493827 ;
weight(1,5) = 0.790123456790123 ; weight(1,6) = 0.493827160493827 ;
weight(1,7) = 0.308641975308642 ; weight(1,8) = 0.493827160493827 ;
weight(1,9) = 0.308641975308642 ;

end
