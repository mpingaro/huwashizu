function [misesnod,misesel] = computemises(stress,element)

st = stress';

misesnod = (st(:,1).^2+st(:,3).^2+3*st(:,2).^2-st(:,1).*st(:,3)).^(1/2);

%sforzi nodali
sxx = stress(1,:)';
sxy = stress(2,:)';
syy = stress(3,:)';
%sforzi elementari
s_xx = 1/3*(sxx(element(:,1))+ sxx(element(:,2))+ sxx(element(:,3)));
s_xy = 1/3*(sxy(element(:,1))+ sxy(element(:,2))+ sxy(element(:,3)));
s_yy = 1/3*(syy(element(:,1))+ syy(element(:,2))+ syy(element(:,3)));  
misesel = s_xx.^2+s_yy.^2-s_xx.*s_yy+3*s_xy.^2; 