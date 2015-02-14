function [jac,inv_jac,djac] = jacobian_tetra(P)

%% Jacobian matrix

jac = [P(2,1)-P(1,1), P(3,1)-P(1,1), P(4,1)-P(1,1);
       P(2,2)-P(1,2), P(3,2)-P(1,2), P(4,2)-P(1,2);
       P(2,3)-P(1,3), P(3,3)-P(1,3), P(4,3)-P(1,3)];

djac = det(jac);
inv_jac = inv(jac)';

end