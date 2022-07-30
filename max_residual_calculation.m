function [max_residual] = max_residual_calculation(Nx,Ny,Nz,u,uOld,v,vOld,w,wOld)

     for k = 1:Nz+1
uX = u{k}(2:Nx,2:Ny+1);
uXO = uOld{k}(2:Nx,2:Ny+1);
vX = v{k}(2:Nx+1,2:Ny);
vXO = vOld{k}(2:Nx+1,2:Ny);
wX = w{k}(2:Nx+1,2:Ny+1);
wXO = wOld{k}(2:Nx+1,2:Ny+1);
cmax1=max(max(abs((uX-uXO))));
cmax2=max(max(abs((vX-vXO))));
cmax3=max(max(abs((wX-wXO))));
cmax = max([cmax1, cmax2, cmax3]);
max_residual_vector(1,k)=cmax;
max_residual = max(max_residual_vector);
     end

end

