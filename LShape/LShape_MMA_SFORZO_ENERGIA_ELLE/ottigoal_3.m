function [f0val,df0dx,fval,dfdx] = ottigoal_3(x)

%function [goal, sensi] = ottigoal_3(x)

global nelem mc2 KELE young youngmin W D ngdls KASSEM spost dx dy VMAX DYOUNG
global strain DELE stress IOPTION IFILTER ISTRESS YOUNG pexp qexp HH HS F FREEDOFS m n FNEW
global NUMFREEDISP SY pois nnod element mc2 coordinates
global NSTRESS ISTRESS STRESSET

YOUNG = youngmin+ x.^pexp*(young-youngmin);
DYOUNG = pexp*x.^(pexp-1)*(young-youngmin);

%constraint and its gradients
fval = zeros(m,1); dfdx = zeros(m,n);
fval(1) = [0.25*dx*dy*sum(x)/VMAX - 1];
dfdx(1,:) = 0.25*dx*dy/VMAX*ones(1,n);

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM] = assembly_goal_11();

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
    %if ISTRESS == 1
    %    stress = reshape(stress',3,[]);
    %    %sforzi nodali
    %    sxx = stress(1,:)';
    %    sxy = stress(2,:)';
    %    syy = stress(3,:)';
    %    %sforzi elementari
    %    s_xx = 1/3*(sxx(element(:,1))+ sxx(element(:,2))+ sxx(element(:,3)));
    %    s_xy = 1/3*(sxy(element(:,1))+ sxy(element(:,2))+ sxy(element(:,3)));
    %    s_yy = 1/3*(syy(element(:,1))+ syy(element(:,2))+ syy(element(:,3)));  
    %    smises2 = s_xx.^2+s_yy.^2-s_xx.*s_yy+3*s_xy.^2; 
    %    
    %    %sensitività sforzi elementari
    %    [usense,sigmasense] = compute_usens(spost,x);
    %    spxx = sigmasense(1:3:end,:);
    %    spxy = sigmasense(2:3:end,:);
    %    spyy = sigmasense(3:3:end,:);    
    %    sprime_xx = 1/3*(spxx(element(:,1),:)+ spxx(element(:,2),:)+ spxx(element(:,3),:));
    %    sprime_xy = 1/3*(spxy(element(:,1),:)+ spxy(element(:,2),:)+ spxy(element(:,3),:));
    %    sprime_yy = 1/3*(spyy(element(:,1),:)+ spyy(element(:,2),:)+ spyy(element(:,3),:));  
    %    smises2 = s_xx.^2+s_yy.^2-s_xx.*s_yy+3*s_xy.^2;         
    %end    
    goal = 0;
    sensi = zeros(nelem,1);
    if IOPTION == 2 
        goal = 0;     
        icountstress = 0;
        for k=1:nelem
            indice = rem(k,4);
            if indice == 0;
                indice = 4;
            end
            epsielement =  strain(mc2(k,:));
            deviaelement = deviatore(epsielement);
            kelement = YOUNG(k)*KELE(:,:,indice);
            goale = satura(kelement*epsielement,epsielement);
            goal = goal + goale;
            sensi(k) = -DYOUNG(k)*goale;
            if ISTRESS == 1 && ismember(k,STRESSET)
                icountstress = icountstress+1;
                [~,epsisense] = compute_epsisense(spost,x);
                epsisenselement =  epsisense(mc2(k,:),:);   
                deviasenselment = deviatore(epsisenselement);        
                fval(1+icountstress) = (YOUNG(k)*pois(4))^2*satura(deviaelement,deviaelement) - (SY*x(k)^qexp)^2/6;
                dfdx(1+icountstress,:) = 2*(YOUNG(k)*pois(4))^2*satura(repmat(deviaelement,1,n),deviasenselment);              
                dfdx(1+icountstress,icountstress) = dfdx(1+icountstress,icountstress) + ...
                    2*pois(4)^2*YOUNG(k)*DYOUNG(k)*satura(deviaelement,deviaelement)- 2/6*qexp*(x(k)^(2*qexp-1))*SY^2;
            end
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
    sensi = HH*(x(:).*sensi(:))./HS./max(1e-3,x(:));
end

f0val = goal;
df0dx = sensi;

%fval = [0.25*dx*dy*sum(x)/VMAX - 1
%    smises2-(x.^(2*qexp))*SY^2];

%s_xx = repmat(s_xx,1,n);
%s_xy = repmat(s_xy,1,n);
%s_yy = repmat(s_yy,1,n);
%   
%dfdx = zeros(m,n);
%dfdx(1,:) = 0.25*dx*dy/VMAX*ones(1,n);
%dfdx(2:end,:) = [2*s_xx.*sprime_xx + 2*s_yy.*sprime_yy + 6*s_xy.*sprime_xy ...
%    - s_xx.*sprime_yy - s_yy.*sprime_xx - 2*qexp*repmat((x.^(2*qexp-1))*SY^2,1,n)]';  
end



 