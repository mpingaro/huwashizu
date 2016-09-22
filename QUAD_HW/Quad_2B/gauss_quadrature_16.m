% by Marco Pingaro & Paolo Venini

function [weight,gauss_x,gauss_y] = gauss_quadrature_16()

%% Gauss Quadrature (16 points)
%
qp = [-0.8611363115940525, ...
    -0.3399810435848563, 0.8611363115940525, 0.3399810435848563];

wp = [0.3478548451374544, ... 
    0.6521451548625460, 0.3478548451374544, 0.6521451548625460];

npg = 16; 
gauss_x = zeros(1,npg);
gauss_y = zeros(1,npg);
weight = zeros(1,npg);
for i = 1:4
    for j = 1:4
        gauss_x(1, 4*(i-1)+j) = qp(i);
        gauss_y(1, 4*(i-1)+j) = qp(j);
        weight(1, 4*(i-1)+j) = wp(i)*wp(j);
    end
end

return