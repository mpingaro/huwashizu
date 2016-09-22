% by Marco Pingaro & Paolo Venini

function [KASSEM,D,W,B,M,K] = assembly_beam(coordinates,element,mc,mc2,lambda,alpha,mu,nelem,ngdlu,ngdls)

%% GLOBAL STIFF MATRIX (PREALLOCATION)
A = sparse(ngdlu,ngdlu);
B = sparse(ngdls,ngdlu);
W = sparse(ngdls,ngdlu);
K = sparse(ngdls,ngdls);
M = sparse(ngdls,ngdls);
D = sparse(ngdls,ngdls);

%% ASSEMBLY GLOBAL MATRIX
for k = 1:nelem
    P(1,[1 2]) = coordinates(element(k,1),[1 2]);
    P(2,[1 2]) = coordinates(element(k,2),[1 2]);
    P(3,[1 2]) = coordinates(element(k,3),[1 2]);
    P(4,[1 2]) = coordinates(element(k,4),[1 2]);
    [AELEM,BELEM,KELEM,MELEM,WELEM,DELEM] = reddy_element(P,lambda,mu);
    % ASSEMBLY A, B, W;
    for i=1:12
        for j=1:12
           A(mc(k,i), mc(k,j))   = A(mc(k,i), mc(k,j))   + AELEM(i,j);
           K(mc2(k,i), mc2(k,j)) = K(mc2(k,i), mc2(k,j)) + KELEM(i,j);
           M(mc2(k,i), mc2(k,j)) = M(mc2(k,i), mc2(k,j)) + MELEM(i,j); 
           B(mc2(k,i), mc(k,j))  = B(mc2(k,i), mc(k,j))  + BELEM(i,j);
           W(mc2(k,i), mc(k,j))  = W(mc2(k,i), mc(k,j))  + WELEM(i,j);
        end
        D(mc2(k,i),mc2(k,i)) = D(mc2(k,i),mc2(k,i)) + DELEM(i,i);
    end
end

D = inv(D);
%% GLOBAL SYSTEM
KASSEM = alpha.*A - alpha.*(B'*D*W + W'*D*B) + W'*D*(K + alpha.*M)*D*W;

end
