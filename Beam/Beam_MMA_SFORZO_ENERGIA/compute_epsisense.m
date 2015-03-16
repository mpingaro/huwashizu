% by Marco Pingaro & Paolo Venini

function [usense,epsisense] = compute_epsisense(spost,x)

global KASSEM FIXEDDOFS FREEDOFS D DYOUNG ngdlu ngdls nelem KELE mc2 W D
global A B M K alpha

LK  = length(KASSEM);
LLK = length(K); 

usense = sparse(LK,nelem);
prodku = sparse(LK,nelem);
%sigmasense = sparse(LLK,nelem);
epsisense = sparse(LLK,nelem);

for k=1:nelem
    indice = rem(k,4); 
    if indice == 0 
        indice = 4;
    end
    KPRIME = sparse(ngdls,ngdls);
    KPRIME(mc2(k,:),mc2(k,:)) = DYOUNG(k)*KELE(:,:,indice);
    %sigmasense(:,k) = (D*KPRIME*D*W)*spost;
    prodku(:,k)= W'*(D*KPRIME*D*W)*spost;     
end

%% SOLVE LINEAR SYSTEMS
usense(FREEDOFS,:) = -KASSEM(FREEDOFS,FREEDOFS)\prodku(FREEDOFS,:);
epsisense = D*W*usense;
%sigmasense = sigmasense + D*(-alpha*B+(K+alpha*M)*D*W)*usense;

end

