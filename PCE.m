function [pPrime] = PCE(Nx,Ny,Nz,dx,dy,dz,dU,dV,dW,uStar,vStar,wStar,pPrime,n,maxNMiter)

 % solve pressure correction equation (PCE)
 % setup coefficients
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny+1
 % boundary coefficients
 aPE{k}(i,j) = dU{k}(i,j)*dy*dz;
 aPW{k}(i,j) = dU{k}(i-1,j)*dy*dz;
 aPN{k}(i,j) = dV{k}(i,j)*dx*dz;
 aPS{k}(i,j) = dV{k}(i,j-1)*dx*dz;
 aPF{k}(i,j) = dW{k}(i,j)*dx*dy;
 aPB{k}(i,j) = dW{k-1}(i,j)*dx*dy;
 % central coefficient
 aPP{k}(i,j) = aPE{k}(i,j)+aPW{k}(i,j)+aPN{k}(i,j)+aPS{k}(i,j)+aPF{k}(i,j)+aPB{k}(i,j);
 % RHS value
 bPP{k}(i,j)=(uStar{k}(i-1,j)-uStar{k}(i,j))*dy*dz+(vStar{k}(i,j-1)-vStar{k}(i,j))*dx*dz...
     + (wStar{k-1}(i,j)-wStar{k}(i,j))*dx*dy;
 end
 end
 end
 
 % fix pressure to zero at bottom left cell
 if n == 1
 i = 2;
 j = 2;
 k = 2;
 aPP{k}(i,j) = 1.0;
 aPE{k}(i,j) = 0.0;
 aPW{k}(i,j) = 0.0;
 aPN{k}(i,j) = 0.0;
 aPS{k}(i,j) = 0.0;
 bPP{k}(i,j) = 0.0;
 end
 
 for m=1:Nz+2
    pPrime{m}=zeros(Nx+2,Ny+2);
 end
 
 for iter = maxNMiter
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny+1
 pPrime{k}(i,j) =(aPE{k}(i,j)*pPrime{k}(i+1,j)+ aPW{k}(i,j)*pPrime{k}(i-1,j)...
     +aPN{k}(i,j)*pPrime{k}(i,j+1)+aPS{k}(i,j)*pPrime{k}(i,j-1)+...
     +aPF{k}(i,j)*pPrime{k+1}(i,j)+aPB{k}(i,j)*pPrime{k-1}(i,j)+bPP{k}(i,j))/aPP{k}(i,j);
 end
 end
 end
 end

end

