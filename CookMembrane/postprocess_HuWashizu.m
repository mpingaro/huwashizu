% by Marco Pingaro & Paolo Venini

function [defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha)

% Deformed mesh
nnod = size(coordinates,1);
Ux = spost(1:2:2*nnod); 
Uy = spost(2:2:2*nnod);
defo =coordinates + [Ux, Uy];

formatSpec = 'Massimo spostamento nella direzione y nel vertice A %12.8f\n';
max_prost = full(Uy(end));
fprintf (formatSpec, max_prost);

% Recovery Strain
strain = D*W*spost;
% Recovery Stress
stress = D*(-alpha*B*spost + (K + alpha*M)*strain);

end