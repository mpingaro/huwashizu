function [JJ,DJ] = jacobian(point)

%[weight,gauss_x,gauss_y] = gauss_quadrature();
[weight, gauss] = GaussQuad2D(4,4);
gauss_x = gauss(:,1);
gauss_y = gauss(:,2);

n_qp = size(weight,1); % prima 2
JJ   = zeros(2,2,n_qp);
DJ   = zeros(1,n_qp);

for i = 1:n_qp
    %x = gauss_x(1,i);
    %y = gauss_y(1,i);
  
    x = gauss_x(i,1);
    y = gauss_y(i,1);
   
    % Shape function gradient
    % deriv along first direction
    grad(1,1:4) = [-(1-y), 1-y, 1+y, -(1+y)].*0.25 ;
    % deriv along second direction
    grad(2,1:4) = [-(1-x), -(1+x), 1+x, 1-x].*0.25 ;
   
    % Jacobian Matrix 
%     J(1,1) = sum( grad(1,1:4)*point(1:4,1) ) ; % x_u
%     J(1,2) = sum( grad(2,1:4)*point(1:4,1) ) ; % x_v
%     J(2,1) = sum( grad(1,1:4)*point(1:4,2) ) ; % y_u
%     J(2,2) = sum( grad(2,1:4)*point(1:4,2) ) ; % y_v
    
    J(1,1) = grad(1,:)*point(:,1) ; % x_u
    J(1,2) = grad(2,:)*point(:,1) ; % x_v
    J(2,1) = grad(1,:)*point(:,2) ; % y_u
    J(2,2) = grad(2,:)*point(:,2) ; % y_v
    
    % Determinant of Jacobian Matrix
    DJ(i) = J(1,1)*J(2,2)-J(1,2)*J(2,1) ;
    % Inverse transpose of Jacobian Matrix
    JJ(1,1,i) =  J(2,2)/DJ(i) ;  
    JJ(1,2,i) = -J(1,2)/DJ(i) ;
    JJ(2,1,i) = -J(2,1)/DJ(i) ; 
    JJ(2,2,i) =  J(1,1)/DJ(i) ;

end

return
