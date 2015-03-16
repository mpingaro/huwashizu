% by Marco Pingaro & Paolo Venini

function [FNEW,spost] = solve_HuWashizu_goal_1()

global FORIGINAL KASSEM FIXEDDOFS FORCEDDOFS FREEDOFS F

F = FORIGINAL;

ALLDOFS = 1:length(KASSEM);
FREEDOFS = setdiff(ALLDOFS,FIXEDDOFS(1,:));
F(FORCEDDOFS(1,:)) =  FORCEDDOFS(2,:);
spost=zeros(length(KASSEM),1);
spost(FIXEDDOFS(1,:),1) = FIXEDDOFS(2,:);

%% SOLVE LINEAR SYSTEM
spost(FREEDOFS,1) = KASSEM(FREEDOFS,FREEDOFS)\F(FREEDOFS,1);

%% RETURN UPDATED LOAD
FNEW = F;

end

