% by Marco Pingaro & Paolo Venini

% Dirichlet or Neumann boundary conditions.
function bc = boundary(ndx,ndy,val,side)

if side == 1
    bc(1,:) = 1:2*(ndx+1);
    bc(2,:) = repmat(val,1,ndx+1); 

elseif side == 2
    bc(1,:) = 2*(ndx+1)*ndy+2*ndx*ndy+1:2*(ndx+1)*(ndy+1)+2*ndx*ndy;
    bc(2,:) = repmat(val,1,ndx+1);

elseif side == 3
    bc = zeros(2,2*(ndy+1));
    for i= 1:ndy+1
        bc(1,2*(i-1)+1) = 2*((ndx+1)*(i-1)+1)-1;
        bc(1,2*(i-1)+2) = 2*((ndx+1)*(i-1)+1);
    end
    bc(2,:) = repmat(val,1,ndy+1);
    
elseif side == 4
    bc = zeros(2,2*(ndy+1));
    for i= 1:ndy+1
        bc(1,2*(i-1)+1) =2*(ndx+1)*i-1;
        bc(1,2*(i-1)+2) =2*(ndx+1)*i;
    end
    bc(2,:) = repmat(val,1,ndy+1);
    
else
    fprintf('Warning: index out of bounds');

end

end
