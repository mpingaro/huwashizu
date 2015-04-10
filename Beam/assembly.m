% by Marco Pingaro & Paolo Venini

function [KASSEM,MASSEM,F] = assembly(coordinates,element,mc,mc2,lambda,alpha,mu,rho,g,nelem,ngdlu,ngdls)

%% GLOBAL STIFF MATRIX (PREALLOCATION)
A = sparse(ngdlu,ngdlu);
B = sparse(ngdls,ngdlu);
W = sparse(ngdls,ngdlu);
K = sparse(ngdls,ngdls);
M = sparse(ngdls,ngdls);
D = sparse(ngdls,ngdls);
F = sparse(ngdlu,1);
MM = sparse(ngdlu,ngdlu);

%% ASSEMBLY GLOBAL MATRIX
for k = 1:nelem
    P(1,[1 2]) = coordinates(element(k,1),[1 2]);
    P(1,[3 4]) = coordinates(element(k,2),[1 2]);
    P(1,[5 6]) = coordinates(element(k,3),[1 2]);
    [AELEM,BELEM,KELEM,MELEM,WELEM,DELEM] = reddy_element(P,lambda,mu);
    MASS = mass_matrix(P,rho); 
    load = body_load(P,g);
    % ASSEMBLY A, B, W;
    for i=1:8
        for j=1:8
           A(mc(k,i),mc(k,j)) = A(mc(k,i),mc(k,j)) + AELEM(i,j);
           MM(mc(k,i),mc(k,j)) = MM(mc(k,i),mc(k,j)) + MASS(i,j);
        end
        for kk = 1:9
            B(mc2(k,kk), mc(k,i)) = B(mc2(k,kk),mc(k,i)) + BELEM(kk,i);
            W(mc2(k,kk), mc(k,i)) = W(mc2(k,kk),mc(k,i)) + WELEM(kk,i);
        end
    end
    % ASSEMBLY K, M, D
    for ii = 1:9
        for jj = 1:9
            K(mc2(k,ii),mc2(k,jj)) = K(mc2(k,ii),mc2(k,jj)) + KELEM(ii,jj);
            M(mc2(k,ii),mc2(k,jj)) = M(mc2(k,ii),mc2(k,jj)) + MELEM(ii,jj);
        end
        D(mc2(k,ii),mc2(k,ii)) = D(mc2(k,ii),mc2(k,ii)) + DELEM(ii,ii);
    end
    % ASSEMBLY BODY LOAD VECTOR
    for z = 1:8
        F(mc(k,z),1) = F(mc(k,z),1) + load(z,1);
    end
end

D = inv(D);
%% GLOBAL SYSTEM
%KASSEM = alpha.*A - alpha.*(B'*D*W + W'*D*B) + W'*D*(K + alpha.*M)*D*W;
KASSEM = []; % matrice globale del sistema 
MASSEM = []; % matrice globale delle masse

end
