% by Marco Pingaro & Paolo Venini

function [KASSEM] = assembly_goal_1()

global YOUNG DYOUNG mc2 alpha nelem ngdls
global KELE 
global A B W M D K

%% GLOBAL STIFF MATRIX (PREALLOCATION)
K = sparse(ngdls,ngdls);
KPRIME = sparse(ngdls,ngdls);

%% ASSEMBLY GLOBAL MATRIX
for k = 1:nelem
    indice = rem(k,4); 
    if indice == 0
        indice = 4;
    end
    KELEM = YOUNG(k)*KELE(:,:,indice);
    for ii = 1:9
        for jj = 1:9
            K(mc2(k,ii),mc2(k,jj)) = K(mc2(k,ii),mc2(k,jj)) + KELEM(ii,jj);
        end
    end
end
%D = inv(D);

%% GLOBAL SYSTEM
KASSEM = alpha.*A - alpha.*(B'*D*W + W'*D*B) + W'*D*(K + alpha.*M)*D*W;

end