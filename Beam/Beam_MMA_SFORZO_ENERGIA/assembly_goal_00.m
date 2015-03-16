% by Marco Pingaro & Paolo Venini

function [A,B,W,M,D,F] = assembly_goal_00()

global coordinates element mc g nelem ngdlu 
global AELE BELE WELE MELE DELE D ik jk iik jjk iiik jjjk

F = sparse(ngdlu,1);
%
AELE1 = AELE(:,:,1);
AELE2 = AELE(:,:,2);
AELE3 = AELE(:,:,3);
AELE4 = AELE(:,:,4);
%
BELE1 = BELE(:,:,1);
BELE2 = BELE(:,:,2);
BELE3 = BELE(:,:,3);
BELE4 = BELE(:,:,4);
%
MELE1 = MELE(:,:,1);
MELE2 = MELE(:,:,2);
MELE3 = MELE(:,:,3);
MELE4 = MELE(:,:,4);
%
WELE1 = WELE(:,:,1);
WELE2 = WELE(:,:,2);
WELE3 = WELE(:,:,3);
WELE4 = WELE(:,:,4);
%
DELE1 = DELE(:,:,1);
DELE2 = DELE(:,:,2);
DELE3 = DELE(:,:,3);
DELE4 = DELE(:,:,4);
%

%% GLOBAL STIFF MATRIX (PREALLOCATION)
soglia = 1.e-11;

%
sk1 = reshape(AELE1(:)*ones(1,nelem/4),64*nelem/4,1);
sk2 = reshape(AELE2(:)*ones(1,nelem/4),64*nelem/4,1);
sk3 = reshape(AELE3(:)*ones(1,nelem/4),64*nelem/4,1);
sk4 = reshape(AELE4(:)*ones(1,nelem/4),64*nelem/4,1);
sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 
A = sparse(iiik,jjjk,sk); A = (A+A')/2;
A(abs(A)<soglia)=0; 
%
sk1 = reshape(MELE1(:)*ones(1,nelem/4),81*nelem/4,1);
sk2 = reshape(MELE2(:)*ones(1,nelem/4),81*nelem/4,1);
sk3 = reshape(MELE3(:)*ones(1,nelem/4),81*nelem/4,1);
sk4 = reshape(MELE4(:)*ones(1,nelem/4),81*nelem/4,1);
clear sk
sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 
M = sparse(ik,jk,sk); M = (M+M')/2;
M(abs(M)<soglia)=0; 
%
sk1 = reshape(DELE1(:)*ones(1,nelem/4),81*nelem/4,1);
sk2 = reshape(DELE2(:)*ones(1,nelem/4),81*nelem/4,1);
sk3 = reshape(DELE3(:)*ones(1,nelem/4),81*nelem/4,1);
sk4 = reshape(DELE4(:)*ones(1,nelem/4),81*nelem/4,1);
clear sk
sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 
D = sparse(diag(diag(sparse(ik,jk,sk)))); %D = (D+D')/2;
D(abs(D)<soglia)=0; 
%
sk1 = reshape(BELE1(:)*ones(1,nelem/4),72*nelem/4,1);
sk2 = reshape(BELE2(:)*ones(1,nelem/4),72*nelem/4,1);
sk3 = reshape(BELE3(:)*ones(1,nelem/4),72*nelem/4,1);
sk4 = reshape(BELE4(:)*ones(1,nelem/4),72*nelem/4,1);
clear sk
sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 
B = sparse(iik,jjk,sk); 
B(abs(B)<soglia)=0; 
%
sk1 = reshape(WELE1(:)*ones(1,nelem/4),72*nelem/4,1);
sk2 = reshape(WELE2(:)*ones(1,nelem/4),72*nelem/4,1);
sk3 = reshape(WELE3(:)*ones(1,nelem/4),72*nelem/4,1);
sk4 = reshape(WELE4(:)*ones(1,nelem/4),72*nelem/4,1);
clear sk
sk(1:length(sk1)) = sk1; 
sk(length(sk1)+1:2*length(sk1)) = sk2; 
sk(2*length(sk1)+1:3*length(sk1)) = sk3; 
sk(3*length(sk1)+1:4*length(sk1)) = sk4; 
W = sparse(iik,jjk,sk); 
W(abs(W)<soglia)=0; 
clear sk

%% ASSEMBLY LOAD VECTOR
for k = 1:nelem
    P(1,[1 2]) = coordinates(element(k,1),[1 2]);
    P(1,[3 4]) = coordinates(element(k,2),[1 2]);
    P(1,[5 6]) = coordinates(element(k,3),[1 2]);
    carico = body_load(P,g);    
    for z = 1:8
        F(mc(k,z),1) = F(mc(k,z),1) + carico(z,1);
    end
end

D = inv(D);

clear AELE BELE WELE MELE DELE
clear AELE1 AELE2 AELE3 AELE4 BELE1 BELE2 BELE3 BELE4
clear WELE1 WELE2 WELE3 WELE4 MELE1 MELE2 MELE3 MELE4
clear DELE1 DELE2 DELE3 DELE4

end