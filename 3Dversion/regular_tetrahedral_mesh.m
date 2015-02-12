function [T,V] = regular_tetrahedral_mesh(lx,ly,lz,nx,ny,nz)
  % REGULAR_TETRAHEDRAL_MESH
  %
  % [T,V] = regular_tetrahedral_mesh(lx,ly,lz,nx,ny,nz)
  %
  % Generates a regular tetrahedral mesh with dimensions (nx,ny,nz)
  %
  % Input:
  %   nx  number of points in x direction on grid
  %   ny  number of points in y direction on grid
  %   nz  number of points in z direction on grid
  % Output:
  %   T  tetrahedra list of indices into V
  %   V  list of vertex coordinates in 3D
  %
  
  % Create a grid
  [x,y,z] = meshgrid(linspace(0,lx,nx),linspace(0,ly,ny),linspace(0,lz,nz));
  V = [x(:) y(:) z(:)];
  % meshgrid flips x and y ordering
  idx = reshape(1:prod([ny,nx,nz]),[ny,nx,nz]);
  v1 = idx(1:end-1,1:end-1,1:end-1);v1=v1(:);
  v2 = idx(1:end-1,2:end,1:end-1);v2=v2(:);
  v3 = idx(2:end,1:end-1,1:end-1);v3=v3(:);
  v4 = idx(2:end,2:end,1:end-1);v4=v4(:);
  v5 = idx(1:end-1,1:end-1,2:end);v5=v5(:);
  v6 = idx(1:end-1,2:end,2:end);v6=v6(:);
  v7 = idx(2:end,1:end-1,2:end);v7=v7(:);
  v8 = idx(2:end,2:end,2:end);v8=v8(:);
  T = [ ...
    v1  v3  v8  v7; ...
    v1  v8  v5  v7; ...
    v1  v3  v4  v8; ...
    v1  v4  v2  v8; ...
    v1  v6  v5  v8; ...
    v1  v2  v6  v8];
end