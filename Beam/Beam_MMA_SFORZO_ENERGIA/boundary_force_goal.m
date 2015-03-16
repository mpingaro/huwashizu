% by Marco Pingaro & Paolo Venini

% Neumann boundary conditions (concentrate force)
function bc = boundary_force_goal(ndx,ndy,val,side)

if side == 1
    bc(1,:) = [1,2];
    bc(2,:) = val; 

elseif side == 2
    bc(1,:) = [2*(ndx+1)-1, 2*(ndx+1)];
    bc(2,:) = val;

elseif side == 3
    nnod = (ndx+1)*(ndy+1)+ndx*ndy-ndx;
    bc(1,:) = [2*nnod-1, 2*nnod];
    bc(2,:) = val;
    
elseif side == 4
    %nnod = (ndx+1)*(ndy+1)+ndx*ndy;
    nnod = (ndx+1)*(ndy+1)-ndx;
    bc(1,:) = [2*nnod-1, 2*nnod];
    bc(2,:) = val;
    
else
    fprintf('Warning: index out of bounds');

end

end
