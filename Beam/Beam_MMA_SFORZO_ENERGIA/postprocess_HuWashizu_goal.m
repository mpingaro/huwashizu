% by Marco Pingaro & Paolo Venini

function [strain,stress] = postprocess_HuWashizu_goal()

global spost D W B M K alpha

% Recovery Strain
strain = D*W*spost;
% Recovery Stress
stress = D*(-alpha*B*spost + (K + alpha*M)*strain);

end