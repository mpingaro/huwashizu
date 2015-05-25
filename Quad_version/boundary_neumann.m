% by Marco Pingaro & Paolo Venini

% Dirichlet or Neumann boundary conditions.
function bc = boundary_neumann(ndx,ndy,val,side)

if side == 1
    val = val/ndx;
    bc(1,:) = 1:2*(ndx+1);
    bc(2,:) = repmat(val,1,ndx+1);
    bc(2,[1, 2]) = val([1 2])/2;
    bc(2,[end-1, end]) = val([1 2])/2; 

elseif side == 2
    val = val/ndx;
    bc(1,:) = 2*(ndx+1)*ndy+1:2*(ndx+1)*(ndy+1);
    bc(2,:) = repmat(val,1,ndx+1);
    bc(2,[1, 2]) = val([1 2])/2;
    bc(2,[end-1, end]) = val([1 2])/2;

elseif side == 3
    val = val/ndy;
    bc = zeros(2,2*(ndy+1));
    for i= 1:ndy+1
        bc(1,2*(i-1)+1) = 2*((ndx+1)*(i-1)+1)-1;
        bc(1,2*(i-1)+2) = 2*((ndx+1)*(i-1)+1);
    end
    bc(2,:) = repmat(val,1,ndy+1);
    bc(2,[1, 2]) = val([1 2])/2;
    bc(2,[end-1, end]) = val([1 2])/2;
    
elseif side == 4
    val = val/ndy;
    bc = zeros(2,2*(ndy+1));
    for i= 1:ndy+1
        bc(1,2*(i-1)+1) =2*(ndx+1)*i-1;
        bc(1,2*(i-1)+2) =2*(ndx+1)*i;
    end
    bc(2,:) = repmat(val,1,ndy+1);
    bc(2,[1, 2]) = val([1 2])/2;
    bc(2,[end-1, end]) = val([1 2])/2;
    
else
    fprintf('Warning: index out of bounds');

end

end