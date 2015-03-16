% by Marco Pingaro & Paolo Venini

function spost = solve_HuWashizu(KASSEM,F,ndx,ndy,bcn,fn,bct,ft,bcd,ud)

%% BOUNDARY IMPOSITION
% Neumann boundary conditions (edges)
for iside=bcn
    ncont = 1;
    bc = boundary_neumann(ndx,ndy,fn(iside,:),iside);
    for gdl=bc(1,:)
        F(gdl,1) = F(gdl,1) + bc(2,ncont);
        ncont = ncont+1;
    end
end
% Neumann boundary conditions (vertex)
for iside=bct
    bt = boundary_force(ndx,ndy, ft(iside,:),iside);
    F(bt(1,1),1) = F(bt(1,1),1) + bt(2,1);
    F(bt(1,2),1) = F(bt(1,2),1) + bt(2,2);
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
