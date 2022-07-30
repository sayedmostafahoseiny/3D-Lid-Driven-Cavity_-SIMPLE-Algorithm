function [Nx,Ny,Nz,dx,x,dy,y,dz,z] = Mesh(L1,L2,L3)
 
% initialize grid matrix (Mesh)
Nx = 25; % number of cells in x direction
dx = L1/Nx;
x = 0:dx:L1;
Ny = 25; % number of cells in y direction
dy = L2/Ny;
y = 0:dy:L2;
Nz = 11;
dz = L3/Nz;
z = 0:dz:L3;

end % End of function

