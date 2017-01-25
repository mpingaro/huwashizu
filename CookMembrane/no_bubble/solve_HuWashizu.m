% by Marco Pingaro & Paolo Venini

function [spost,defo,strain,stress] = solve_HuWashizu(coordinates,ndx,ndy,f_traction,KASSEM,ngdlu,D,W,B,M,K,alpha)

%% BOUNDARY IMPOSITION 
%     Trave incastrata e trazione sul lato verticale destro

% Calcoli indici bordo incastrato e caricato
tra = zeros(ndy+1,1);
inc = zeros(2*(ndy+1),1);
for i= 1:ndy+1
    tra(i,1) = 2*(ndx+1)*i;
    inc(2*(i-1)+1,1) = 2*((ndx+1)*(i-1)+1)-1;
    inc(2*(i-1)+2,1) = 2*((ndx+1)*(i-1)+1);
end

% Imposizione vincoli
for I = 1:size(inc,1)
    KASSEM(inc(I,1), :) = 0;
    KASSEM(inc(I,1),inc(I,1)) = 1;
end

f_traction = f_traction/(2*ndy);
load = sparse(ngdlu,1);
for J = 1:size(tra,1)
    if J==1 || J == size(tra,1)
        load(tra(J,1),1) = f_traction;
    else
        load(tra(J,1),1) = 2*f_traction;
    end
end

% Solution
spost = KASSEM\load;
[defo,strain,stress] = postprocess_HuWashizu(coordinates,spost,D,W,B,M,K,alpha);

end
