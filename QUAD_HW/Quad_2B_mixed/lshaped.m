% Created by Marco Pingaro & Paolo Venini

function [coordinates,element,nnod,nelem] = lshaped(length,heigth,ndx,ndy)

%% Coordinates
%%
px = length/ndx;
py = heigth/ndy;
npoint = (2*ndx+1)*(ndy+1);
xcor = 0:px:2*length;
ycor = 0:py:heigth;

y = repmat(ycor,2*ndx+1,1);
x = repmat(xcor,ndy+1,1);

coord_b = zeros(npoint,2);
coord_b(:,1) = reshape(x',npoint,1);
coord_b(:,2) = reshape(y,npoint,1);
%%
npoint = (ndx+1)*ndy;
xcor = length:px:2*length;
ycor = heigth+py:py:2*heigth;

y = repmat(ycor, ndx+1,1);
x = repmat(xcor, ndy,1);

coord_u = zeros(npoint,2);
coord_u(:,1) = reshape(x',npoint,1);
coord_u(:,2) = reshape(y,npoint,1);

coordinates = [coord_b; coord_u];
nnod = size(coordinates,1);

%% Element
nel_b = 2*ndx*ndy;
element_b = zeros(nel_b,4);
for i =1:ndy
    for j = 1:2*ndx
        % Element
        element_b(2*ndx*(i-1)+j,1) = (2*ndx+1)*(i-1)+j;
        element_b(2*ndx*(i-1)+j,2) = element_b(2*ndx*(i-1)+j,1)+1;
        element_b(2*ndx*(i-1)+j,3) = (2*ndx+1)*i+j+1;
        element_b(2*ndx*(i-1)+j,4) = element_b(2*ndx*(i-1)+j,3)-1;
    end
end

nel_u = ndx*ndy;
element_u = zeros(nel_u,4);
for i =1:ndy
    for j = 1:ndx
        % Element
        element_u(ndx*(i-1)+j,1) = (ndx+1)*(i-1)+(2*ndx+1)*ndy+ndx+j;
        element_u(ndx*(i-1)+j,2) = element_u(ndx*(i-1)+j,1)+1;
        element_u(ndx*(i-1)+j,3) = (ndx+1)*i+(2*ndx+1)*ndy+ndx+j+1;
        element_u(ndx*(i-1)+j,4) = element_u(ndx*(i-1)+j,3)-1;
    end
end

element = [element_b;element_u];
nelem = nel_b+nel_u;

end