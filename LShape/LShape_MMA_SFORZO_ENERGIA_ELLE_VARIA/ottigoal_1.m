function [goal, sensi] = ottigoal_1(x)

global nnod nelem mc2 KELE young youngmin W D ngdls KASSEM spost

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM] = assembly_goal_1(x);

%% SOLVE
[FNEW,spost] = solve_HuWashizu_goal_1();

goal = abs(full(FNEW'*spost));

%[[1:length(FNEW)]' full(FNEW) full(spost)]
%size(FNEW), size(spost)

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


%% POST PROCESSING
%[defo,strain,stress] = postprocess_HuWashizu_goal(coordinates,spost,D,W,B,M,K,alpha);

end

