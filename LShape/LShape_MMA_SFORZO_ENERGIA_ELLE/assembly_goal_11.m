% by Marco Pingaro & Paolo Venini

function [KASSEM] = assembly_goal_11()

global YOUNG alpha nelem 
global A B W M D K ik jk KELE1 KELE2 KELE3 KELE4 KASSEM1

YOUNGSPARSE1 = YOUNG(1:4:end);
YOUNGSPARSE2 = YOUNG(2:4:end);
YOUNGSPARSE3 = YOUNG(3:4:end);
YOUNGSPARSE4 = YOUNG(4:4:end);

sk1 = reshape(KELE1(:)*YOUNGSPARSE1(:)',81*nelem/4,1);
sk2 = reshape(KELE2(:)*YOUNGSPARSE2(:)',81*nelem/4,1);
sk3 = reshape(KELE3(:)*YOUNGSPARSE3(:)',81*nelem/4,1);
sk4 = reshape(KELE4(:)*YOUNGSPARSE4(:)',81*nelem/4,1);

sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 

K = sparse(ik,jk,sk); K = (K+K')/2;

clear sk1 sk2 sk3 sk4

%% GLOBAL SYSTEM
%KASSEM = alpha.*A - alpha.*(B'*D*W + W'*D*B) + W'*D*(K + alpha.*M)*D*W;
KASSEM = KASSEM1 + W'*D*K*D*W;

end