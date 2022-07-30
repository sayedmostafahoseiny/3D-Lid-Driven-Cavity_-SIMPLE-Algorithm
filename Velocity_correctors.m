function [uPrime,vPrime,wPrime] = Velocity_correctors(Nx,Ny,Nz,uPrime,dU,pPrime,vPrime,dV,wPrime,dW)

% calculate velocity corrections
% u corrections
 for k = 2:Nz+1
 for i = 2:Nx
 for j = 2:Ny+1
 uPrime{k}(i,j)=dU{k}(i,j)*(pPrime{k}(i,j)-pPrime{k}(i+1,j));
 end
 end
 end
 % v corrections
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny
 vPrime{k}(i,j)=dV{k}(i,j)*(pPrime{k}(i,j)-pPrime{k}(i,j+1));
 end
 end
 end
 % w corrections
 for k = 2:Nz
 for i = 2:Nx+1
 for j = 2:Ny+1
 wPrime{k}(i,j)=dW{k}(i,j)*(pPrime{k}(i,j)-pPrime{k+1}(i,j));
 end
 end
 end
 
end

