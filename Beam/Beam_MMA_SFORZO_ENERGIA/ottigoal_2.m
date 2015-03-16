function [f0val,df0dx,fval,dfdx] = ottigoal_2(x)

%function [goal, sensi] = ottigoal_2(x)

global nelem mc2 KELE young youngmin W D ngdls KASSEM spost dx dy VMAX
global strain DELE stress IOPTION IFILTER YOUNG pexp HH F FREEDOFS m n FNEW
global NUMFREEDISP

YOUNG = youngmin+ x.^pexp*(young-youngmin);
DYOUNG = pexp*x.^(pexp-1)*(young-youngmin);

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM] = assembly_goal_1();

%% SOLVE
[FNEW,spost] = solve_HuWashizu_goal_1();

%% COMPUTE OBJECTIVE FUNCTION
if IOPTION == 1
    goal = abs(full(FNEW'*spost));
    sensi = zeros(nelem,1);
    for k=1:nelem
        indice = rem(k,4); 
        if indice == 0 
            indice = 4;
        end
        KPRIME = sparse(ngdls,ngdls);
        KPRIME(mc2(k,:),mc2(k,:)) = KELE(:,:,indice);
        KPRIME = W'*D*KPRIME*D*W;
        sensi(k) = -FNEW'*(KASSEM\(3*(x(k)^2)*(young-youngmin)*KPRIME*spost));
    end
else
    [strain,stress] = postprocess_HuWashizu_goal();
    goal = 0;
    sensi = zeros(nelem,1);
    if IOPTION == 2 
        goal = 0;
        for k=1:nelem
            indice = rem(k,4);
            if indice == 0;
                indice = 4;
            end
            epsielement =  strain(mc2(k,:));
            %kelement = YOUNG(k)*KELE(:,:,indice)+alpha*MELE(:,:,indice);
            kelement = YOUNG(k)*KELE(:,:,indice);
            goale = 10^3*(kelement*epsielement)'*epsielement;
            goal = goal + goale;
            sensi(k) = -DYOUNG(k)*goale;
        end      
    else
        for k=1:nelem
            indice = rem(k,4); 
            if indice == 0
                indice = 4;
            end
            epsielement =  strain(mc2(k,:));
            stresselement =  stress(mc2(k,:));
            goal = goal + (DELE(:,:,indice)*stresselement)'*epsielement;
        end
        goal = abs(goal);
    end
end

%SENSITIVITY FILTERING
if IFILTER == 1
    sensif = zeros(nelem,1);
    for k=1:nelem
        for h=1:nelem
            sensif(k) = sensif(k) + HH(k,h)*x(h)*sensi(h);
        end
        sensif(k) = sensif(k)/(x(k)*sum(HH(k,:)));
    end
    sensi = sensif;
end

f0val = goal;
df0dx = sensi;
fval = [0.25*dx*dy*sum(x)-VMAX
        KASSEM(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1)-F(FREEDOFS,1)
       -KASSEM(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1)+F(FREEDOFS,1)];
 dfdx = zeros(m,n);
 dfdx(1,:) = 0.25*dx*dy*ones(1,n);

 for k=1:nelem
     indice = rem(k,4);
     if indice == 0 
         indice = 4;
     end
     KPRIME = sparse(ngdls,ngdls);
     KPRIME(mc2(k,:),mc2(k,:)) = KELE(:,:,indice);
     KPRIME = W'*D*KPRIME*D*W;
     %size(3*(x(k)^2)*(young-youngmin)*KPRIME(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1))
     %size(dfdx(2:NUMFREEDISP+1,k))
     %size([3*(x(k)^2)*(young-youngmin)*KPRIME(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1)])
     %pause
     sens = -KASSEM(FREEDOFS,FREEDOFS)\(DYOUNG(k)*KPRIME(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1));
     dfdx(2:NUMFREEDISP+1,k) = [DYOUNG(k)*KPRIME(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1)+ ... 
         KASSEM(FREEDOFS,FREEDOFS)*sens];
     dfdx(NUMFREEDISP+2:2*NUMFREEDISP+1,k) = -[DYOUNG(k)*KPRIME(FREEDOFS,FREEDOFS)*spost(FREEDOFS,1)+ ... 
         KASSEM(FREEDOFS,FREEDOFS)*sens];
 end

end

