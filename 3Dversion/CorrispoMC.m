% by Marco Pingaro & Paolo Venini

function [mc,ngdlu]=CorrispoMC(element,nelem,nnod)

mc = zeros(nelem,15);         % prealloco matrice di corrispondenza
ngdlu = 3*(nelem+nnod);      % number of degree of freedom displacement

mc(:,1) = 3*element(:,1)-2;
mc(:,2) = 3*element(:,1)-1;
mc(:,3) = 3*element(:,1);
mc(:,4) = 3*element(:,2)-2;
mc(:,5) = 3*element(:,2)-1;
mc(:,6) = 3*element(:,2);
mc(:,7) = 3*element(:,3)-2;
mc(:,8) = 3*element(:,3)-1;
mc(:,9) = 3*element(:,3);
mc(:,10) = 3*element(:,4)-2;
mc(:,11) = 3*element(:,4)-1;
mc(:,12) = 3*element(:,4);
mc(:,13) = 3*nnod+1:3:ngdlu-2;
mc(:,14) = 3*nnod+2:3:ngdlu-1;
mc(:,15) = 3*nnod+2:3:ngdlu;

return