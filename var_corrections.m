function [u,v,w,p] = var_corrections(Nx,Ny,Nz,alphaP,alphaU,p,pStar,...
    pPrime,u,uStar,uPrime,v,vStar,vPrime,w,wStar,wPrime)

% p corrections with under-relaxation
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny+1
 p{k}(i,j)= pStar{k}(i,j) + alphaP * pPrime{k}(i,j);
 end
 end
 end
 % u corrections
 for k = 2:Nz+1
 for i = 2:Nx
 for j = 2:Ny+1
 u{k}(i,j)=uStar{k}(i,j)+ alphaU * uPrime{k}(i,j);
 end
 end
 end
 % v corrections
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny
 v{k}(i,j)=vStar{k}(i,j)+ alphaU * vPrime{k}(i,j);
 end
 end
 end
 % w corrections
 for k = 2:Nz
 for i = 2:Nx+1
 for j = 2:Ny+1
 w{k}(i,j)=wStar{k}(i,j)+ alphaU * wPrime{k}(i,j);
 end
 end
 end
 
end

