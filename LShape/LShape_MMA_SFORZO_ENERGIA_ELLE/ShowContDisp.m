% by Marco Pingaro & Paolo Venini

function ShowContDisp(element,coordinates,r)
hold on
for j=1:size(element,1)
    u(1,[1 2 3]) = r(1,element(j,[1 2 3]));
    trisurf([1 2 3],coordinates(element(j,:),1),...
        coordinates(element(j,:),2),...
        u,'facecolor','interp');
end
