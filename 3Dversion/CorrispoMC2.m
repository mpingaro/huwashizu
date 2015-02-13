% by Marco Pingaro & Paolo Venini

function mc2=CorrispoMC2(element,nelem)

mc2 = zeros(nelem,24);         % prealloco matrice di corrispondenza
 
mc2(:,1)  = 6*element(:,1)-5;
mc2(:,2)  = 6*element(:,1)-4;
mc2(:,3)  = 6*element(:,1)-3;
mc2(:,4)  = 6*element(:,1)-2;
mc2(:,5)  = 6*element(:,1)-1;
mc2(:,6)  = 6*element(:,1);

mc2(:,7)  = 6*element(:,2)-5;
mc2(:,8)  = 6*element(:,2)-4;
mc2(:,9)  = 6*element(:,2)-3;
mc2(:,10) = 6*element(:,2)-2;
mc2(:,11) = 6*element(:,2)-1;
mc2(:,12) = 6*element(:,2);

mc2(:,13) = 6*element(:,3)-5;
mc2(:,14) = 6*element(:,3)-4;
mc2(:,15) = 6*element(:,3)-3;
mc2(:,16) = 6*element(:,3)-2;
mc2(:,17) = 6*element(:,3)-1;
mc2(:,18) = 6*element(:,3);

mc2(:,19) = 6*element(:,4)-5;
mc2(:,20) = 6*element(:,4)-4;
mc2(:,21) = 6*element(:,4)-3;
mc2(:,22) = 6*element(:,4)-2;
mc2(:,23) = 6*element(:,4)-1;
mc2(:,24) = 6*element(:,4);

end