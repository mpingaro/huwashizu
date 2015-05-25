% by Marco Pingaro & Paolo Venini

function [mc,ngdlu]=CorrispoMC(element,nelem,nnod)

mc = zeros(nelem,11);         % prealloco matrice di corrispondenza
ngdlu = 2*nnod + 3*nelem;     % number of degree of freedom displacement

mc(:,1)  = 2*element(:,1)-1;
mc(:,2)  = 2*element(:,1);
mc(:,3)  = 2*element(:,2)-1;
mc(:,4)  = 2*element(:,2);
mc(:,5)  = 2*element(:,3)-1;
mc(:,6)  = 2*element(:,3);
mc(:,7)  = 2*element(:,4)-1;
mc(:,8)  = 2*element(:,4);
mc(:,9)  = 2*nnod+1:3:ngdlu-2;
mc(:,10) = 2*nnod+2:3:ngdlu-1;
mc(:,11) = 2*nnod+3:3:ngdlu;

return