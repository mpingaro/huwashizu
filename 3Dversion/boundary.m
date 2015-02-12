% by Marco Pingaro & Paolo Venini

% Dirichlet or Neumann boundary conditions.
function bc = boundary(coordinates,lx,ly,lz,val,side)

if side == 1
    b =find(coordinates(:,1)==0); % Face 1 x ==0
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1)); 
  
elseif side == 2
    b =find(coordinates(:,1)==lx); % Face 2 x == lx
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1));

elseif side == 3
    b =find(coordinates(:,3)==0); % Face 3 z == 0
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1));
    
elseif side == 4
    b =find(coordinates(:,3)==lz); % Face 4 z == lz
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1));
    
elseif side == 5
    b =find(coordinates(:,2)==0); % Face 5 y == 0
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1));
    
elseif side == 6
    b =find(coordinates(:,2)==ly); % Face 6 y == ly
    bb = [3.*b-2, 3.*b-1, 3.*b]; 
    bc(1,:) = reshape(bb',1,[]);
    bc(2,:) = repmat(val,1, size(b,1));
    
else
    fprintf('Warning: index out of bounds');

end

end
