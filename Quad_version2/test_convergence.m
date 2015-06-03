% By Marco Pingaro & Paolo Venini
clear all;
close all;
clc;
%% TEST CONVERGENCE
nx = [2, 4, 8, 16, 32, 64, 128];
ny = [1, 2, 4, 8, 16, 32, 64];

sp = zeros(size(nx,2),1);
nel = zeros(size(nx,2),1);
for i = 1:size(nx,2)
    px = nx(i);
    py = ny(i);
    %nel(i) = px*py;
    nel(i) = py;
    sp(i) = huwashizu_cook(px,py);
end
figure, plot(nel,sp,'-o');