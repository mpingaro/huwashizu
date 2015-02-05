% by Marco Pingaro & Paolo Venini

function spost = solve_HuWashizu(KASSEM,F,ndx,ndy,bcn,fn,bcd,ud)

%% BOUNDARY IMPOSITION
% Neumann boundary conditions
for iside=bcn
    ncont = 1;
    bc = boundary_neumann(ndx,ndy,fn(iside,:),iside);
    for gdl=bc(1,:)
        F(gdl,1) = F(gdl,1) + bc(2,ncont);
        ncont = ncont+1;
    end
end

% Dirichlet boundary conditions
for iside = bcd
    ncont = 1;
    bc = boundary(ndx,ndy,ud(iside,:),iside);
    for gdl=bc(1,:)
        KASSEM(gdl, :) = 0;
        KASSEM(gdl,gdl) = 1;
        F(gdl,1) = bc(2,ncont);
        ncont = ncont + 1;
    end
end

%% SOLVE LINEAR SYSTEM
spost = KASSEM\F;

end