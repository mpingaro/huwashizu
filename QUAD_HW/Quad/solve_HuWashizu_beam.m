% by Marco Pingaro & Paolo Venini

function spost = solve_HuWashizu_beam(KASSEM,coor,h,f,ndx,ndy,ngdlu)

%% BOUNDARY IMPOSITION
F = sparse(ngdlu,1);
dl = h/ndy;
% Neumann boundary conditions (vertex)
bc = zeros(2,2*(ndy+1));
for i= 1:ndy+1
    bc(1,2*(i-1)+1) =2*(ndx+1)*i-1;
    bc(1,2*(i-1)+2) =2*(ndx+1)*i;

    if (i==1)
        fx = -2*f/h * dl/2 + f;
        fy = 0.0;
        bc(2,2*(i-1)+1) = fx*dl;
        bc(2,2*(i-1)+2) = fy*dl;
    elseif (i==ndy+1)
        fx = -2*f/h *(h-dl/2) + f;
        fy = 0.0;
        bc(2,2*(i-1)+1) = fx*dl;
        bc(2,2*(i-1)+2) = fy*dl;
    else
        fx = -2*f/h * coor(bc(1,2*(i-1)+2)/2,2) + f;
        fy = 0.0;
        bc(2,2*(i-1)+1) = fx*dl;
        bc(2,2*(i-1)+2) = fy*dl;
    end
   
end
ncont = 1;
for gdl=bc(1,:)
    F(gdl,1) = bc(2,ncont);
    ncont = ncont + 1;
end

% Dirichlet boundary conditions
bb = zeros(1,ndy);
for i= 2:ndy+1
    bb(1,i-1) = 2*((ndx+1)*(i-1)+1)-1;
end
bd = [1,2,bb];
ncont = 1;
for gdl=bd(1,:)
    KASSEM(gdl, :) = 0.0;
    KASSEM(:, gdl) = 0.0;
    KASSEM(gdl,gdl) = 1.0;
    F(gdl,1) = 0.0;
    ncont = ncont + 1;
end
%% SOLVE LINEAR SYSTEM
spost = KASSEM\F;

end