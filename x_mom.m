function [auP,dU,uStar] = x_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,uStar,dU)

 % STEP 1a: solve x-momentum as uStar
 % setup coefficients
 for k = 2:Nz+1
 for i = 2:Nx
 for j = 2:Ny+1
 % convective flux using central differencing
 uCe = dy*dz*0.5*(uOld{k}(i,j)+uOld{k}(i+1,j));
 uCw = dy*dz*0.5*(uOld{k}(i,j)+uOld{k}(i-1,j));
 uCn = dx*dz*0.5*(vOld{k}(i,j)+vOld{k}(i+1,j));
 uCs = dx*dz*0.5*(vOld{k}(i,j-1)+vOld{k}(i+1,j-1));
 uCf = dx*dy*0.5*(wOld{k}(i,j)+wOld{k}(i+1,j));
 uCb = dx*dy*0.5*(wOld{k-1}(i,j)+wOld{k-1}(i+1,j));
 uDx = (dy*dz)/(dx*Re);
 uDy = (dx*dz)/(dy*Re);
 uDz = (dx*dy)/(dz*Re);
 
 % boundary coefficients using hybrid schemes
 auE_h = [-uCe,(uDx-0.5*uCe),0];
 auW_h = [uCw,(uDx+0.5*uCw),0];
 auN_h = [-uCn,(uDy-0.5*uCn),0];
 auS_h = [uCs,(uDy+0.5*uCs),0];
 auF_h = [-uCf,(uDz-0.5*uCf),0];
 auB_h = [uCb,(uDz+0.5*uCb),0];
 auE{k}(i,j) = max(auE_h);
 auW{k}(i,j) = max(auW_h);
 auN{k}(i,j) = max(auN_h);
 auS{k}(i,j) = max(auS_h);
 auF{k}(i,j) = max(auF_h);
 auB{k}(i,j) = max(auB_h);
 
 % central coefficient
 auP{k}(i,j) = (auE{k}(i,j)+auW{k}(i,j)+auN{k}(i,j)+auS{k}(i,j)+auF{k}(i,j)+auB{k}(i,j));
 
 % set u velocity component of PCE
 dU{k}(i,j) = (dy*dz*alphaU)/auP{k}(i,j);
 end
 end
 end
 % use previous calculation as initial guess in numerical method
 uStar = uOld;
 
 for iter = maxNMiter
 for k = 2:Nz+1
 for i = 2:Nx
 for j = 2:Ny+1
 uStar{k}(i,j) = (alphaU/auP{k}(i,j)) * (auE{k}(i,j)*uStar{k}(i+1,j)...
 + auW{k}(i,j)*uStar{k}(i-1,j) +auN{k}(i,j)*uStar{k}(i,j+1)...
 + auS{k}(i,j)*uStar{k}(i,j-1)+ auF{k}(i,j)*uStar{k+1}(i,j)...
 + auB{k}(i,j)*uStar{k-1}(i,j) - dy*dz*(pStar{k}(i+1,j)-pStar{k}(i,j)))...
 + (1-alphaU)*uOld{k}(i,j);
 end
 end
 end
 end
 
end % End of function 

