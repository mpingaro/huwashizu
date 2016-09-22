% by Marco Pingaro

function coordinates = beam_distorted(length,heigth,ndx,ndy)

%% Coorinates
px = length/ndx;
py = heigth/ndy;
npoint = (ndx+1)*(ndy+1);
xcor = [0:px:length];
ycor = [0:py:heigth];
for i=1:ndx+1
    y(i,:) = ycor;
end
for j=1:ndy+1
    x(:,j) = xcor; 
end
coordinates = zeros(npoint,1);
coordinates(:,1) = reshape(x,npoint,1);
coordinates(:,2) = reshape(y,npoint,1);
% Insert distorsion
dst = zeros(npoint,1);
for j = 1:ndy/2
    for i = 1:ndx+1
        dst(ndx+1+2*(j-1)*(ndx+1)+i,1) = (-1)^(i)*(py/2); %prima /2
    end
end
coordinates(:,2) = coordinates(:,2)+dst;

end
