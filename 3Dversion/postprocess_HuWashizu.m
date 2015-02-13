% by Marco Pingaro & Paolo Venini

function [defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha)

% Deformed mesh
nnod = size(coordinates,1);
Ux = spost(1:3:3*nnod-2); 
Uy = spost(2:3:3*nnod-1);
Uz = spost(3:3:3*nnod);
defo = coordinates + [Ux, Uy, Uz];

% Recovery Strain
strain = D*W*spost;
% Recovery Stress
stress = D*(-alpha*B*spost + (K + alpha*M)*strain);
end