% by Marco Pingaro & Paolo Venini

% Neumann boundary conditions (concentrate force)
function bc = boundary_force(coor,lx,ly,lz,val,side)

if side == 1 % (0,0,0)
    b = find(coor(:,1) == 0 & coor(:,2) == 0 & coor(:,3) == 0);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val; 

elseif side == 2 % (lx,0,0)
    b = find(coor(:,1) == lx & coor(:,2) == 0 & coor(:,3) == 0);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;

elseif side == 3 % (0,ly,0)
    b = find(coor(:,1) == 0 & coor(:,2) == ly & coor(:,3) == 0);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
elseif side == 4 % (lx,ly,0)
    b = find(coor(:,1) == lx & coor(:,2) == ly & coor(:,3) == 0);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
elseif side == 5 % (0,0,lz)
    b = find(coor(:,1) == 0 & coor(:,2) == 0 & coor(:,3) == lz);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
elseif side == 6 % (lx,0,lz)
    b = find(coor(:,1) == lx & coor(:,2) == 0 & coor(:,3) == lz);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
elseif side == 7 % (0,ly,lz)
    b = find(coor(:,1) == 0 & coor(:,2) == ly & coor(:,3) == lz);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
elseif side == 8 % (lx,ly,lz)
    b = find(coor(:,1) == lx & coor(:,2) == ly & coor(:,3) == lz);
    bc(1,:) = [3*b-2,3*b-1,3*b];
    bc(2,:) = val;
    
else
    fprintf('Warning: index out of bounds');

end

end
