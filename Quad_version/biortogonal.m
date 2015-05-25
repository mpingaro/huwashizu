clear all
close all
clc

% syms x y a b c d
% 
% 
% phi(1) = 0.25.*(1-x).*(1-y);
% phi(2) = 0.25.*(1+x).*(1-y);
% phi(3) = 0.25.*(1+x).*(1+y);
% phi(4) = 0.25.*(1-x).*(1+y);
% 
% t = a + b.*x + c.*y + d.*x*y;
% 
% fun1 = phi(1)*t;
% fun2 = phi(2)*t;
% fun3 = phi(3)*t;
% fun4 = phi(4)*t;
% 
% A1 = int(int(fun1,x,-1,1),y,-1,1)
% A2 = int(int(fun2,x,-1,1),y,-1,1)
% A3 = int(int(fun3,x,-1,1),y,-1,1)
% A4 = int(int(fun4,x,-1,1),y,-1,1)


% a - b/3 - c/3 + d/9
% a + b/3 - c/3 - d/9
% a + b/3 + c/3 + d/9
% a - b/3 + c/3 - d/9


K = [1 -1/3 -1/3 1/9;
     1 1/3 -1/3 -1/9;
     1 1/3 1/3 1/9;
     1 -1/3 1/3 -1/9];
 
b = [0; 0; 0; 4];

x = K\b;




