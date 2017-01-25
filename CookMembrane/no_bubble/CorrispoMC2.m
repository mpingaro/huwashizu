% by Marco Pingaro & Paolo Venini

function mc2=CorrispoMC2(element,nelem)

mc2 = zeros(nelem,9);         % prealloco matrice di corrispondenza
 
mc2(:,1) = 3*element(:,1)-2;
mc2(:,2) = 3*element(:,1)-1;
mc2(:,3) = 3*element(:,1);
mc2(:,4) = 3*element(:,2)-2;
mc2(:,5) = 3*element(:,2)-1;
mc2(:,6) = 3*element(:,2);
mc2(:,7) = 3*element(:,3)-2;
mc2(:,8) = 3*element(:,3)-1;
mc2(:,9) = 3*element(:,3);

return
