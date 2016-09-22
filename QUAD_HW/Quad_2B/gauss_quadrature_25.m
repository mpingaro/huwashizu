% by Marco Pingaro

function [weight,gauss_x,gauss_y] = gauss_quadrature_25()

%% Gauss Quadrature (25 points)
qp = [0.0000000000000000, -0.5384693101056831, 0.5384693101056831,...
    -0.9061798459386640, 0.9061798459386640];
wp = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665,...
    0.2369268850561891, 0.2369268850561891];

npg = 25; 
gauss_x = zeros(1,npg);
gauss_y = zeros(1,npg);
for i = 1:5
    for j = 1:5
        gauss_x(1, 5*(i-1)+j) = qp(i);
        gauss_y(1, 5*(i-1)+j) = qp(j);
        weight(1, 5*(i-1)+j) = wp(i)*wp(j);
    end
end

return
