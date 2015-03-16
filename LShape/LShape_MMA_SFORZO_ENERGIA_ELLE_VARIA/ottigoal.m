function goal = ottigoal_1(x)

global nnod

%% ASSEMBLY GLOBAL MATRIX AND GLOBAL STIFFNESS MATRIX
[KASSEM,K] = assembly_goal_1(x);

%% SOLVE
[FNEW,spost] = solve_HuWashizu_goal(KASSEM,F);
Ux = spost(1:2:2*nnod); 
Uy = spost(2:2:2*nnod);

goal = abs(full(FNEW'*spost));

%% POST PROCESSING
%[defo,strain,stress] = postprocess_HuWashizu_goal(coordinates,spost,D,W,B,M,K,alpha);

end

