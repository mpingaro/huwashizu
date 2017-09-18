% By Marco Pingaro & Paolo Venini
clear;
close all;
clc;
%% TEST CONVERGENCE
%nx = [1, 2, 8, 16, 32, 64];
%ny = [1, 2, 8, 16, 32, 64];
nx = 32;
ny = 32;

name = 'cook_1B_type_2_3mu.txt';
f = fopen( name, 'w' );
fprintf(f, 'number of element per side v.s. vertical diaplacement of point A\n' );

sp = zeros(size(nx,2),1);
nel = zeros(size(nx,2),1);
for i = 1:size(nx,2)
    px = nx(i);
    py = ny(i);
    nel(i) = px*py;
    sp(i) = huwashizu_cook(px,py);
    fprintf(f, '%3.0f \t %5.5e \n', py, sp(i)); 
end
fclose(f);
%figure, plot(nel,sp,'-o');
