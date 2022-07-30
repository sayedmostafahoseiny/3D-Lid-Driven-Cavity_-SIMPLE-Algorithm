function [awP,dW,wStar] = z_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,wStar,dW)
    % STEP 1c: solve z-momentum as wStar
 % setup coefficients
 for k = 2:Nz
 for i = 2:Nx+1
 for j = 2:Ny+1
 % convective flux using central differencing
 wCe = dy*dz*0.5*(uOld{k}(i,j)+uOld{k}(i,j+1));
 wCw = dy*dz*0.5*(uOld{k}(i-1,j)+uOld{k}(i-1,j+1));
 wCn = dx*dz*0.5*(vOld{k}(i,j)+vOld{k}(i+1,j));
 wCs = dx*dz*0.5*(vOld{k}(i,j-1)+vOld{k}(i+1,j-1));
 wCf = dx*dy*0.5*(wOld{k}(i,j)+wOld{k+1}(i,j));
 wCb = dx*dy*0.5*(wOld{k}(i,j)+wOld{k-1}(i,j));
 wDx = (dy*dz)/(dx*Re);
 wDy = (dx*dz)/(dy*Re);
 wDz = (dx*dy)/(dz*Re);
 
 % boundary coefficients using hybrid schemes
 awE_h = [-wCe,(wDx-0.5*wCe),0];
 awW_h = [wCw,(wDx+0.5*wCw),0];
 awN_h = [-wCn,(wDy-0.5*wCn),0];
 awS_h = [wCs,(wDy+0.5*wCs),0];
 awF_h = [-wCf,(wDz-0.5*wCf),0];
 awB_h = [wCb,(wDz+0.5*wCb),0];
 awE{k}(i,j) = max(awE_h);
 awW{k}(i,j) = max(awW_h);
 awN{k}(i,j) = max(awN_h);
 awS{k}(i,j) = max(awS_h);
 awF{k}(i,j) = max(awF_h);
 awB{k}(i,j) = max(awB_h);
 
 % central coefficient using applied under-relaxation
 awP{k}(i,j) = (awE{k}(i,j)+awW{k}(i,j)+awN{k}(i,j)+awS{k}(i,j)+awF{k}(i,j)+awB{k}(i,j));
 
 % set u velocity component of PCE
 dW{k}(i,j) = (dx*dy)/awP{k}(i,j);
 end
 end
 end
 % use previous calculation as initial guess in numerical method
 wStar = wOld;
 
 for iter = maxNMiter
 for k = 2:Nz
 for i = 2:Nx+1
 for j = 2:Ny+1
 wStar{k}(i,j) = (alphaU/awP{k}(i,j)) * (awE{k}(i,j)*wStar{k}(i+1,j)...
 + awW{k}(i,j)*wStar{k}(i-1,j) +awN{k}(i,j)*wStar{k}(i,j+1)...
 + awS{k}(i,j)*wStar{k}(i,j-1)+ awF{k}(i,j)*wStar{k+1}(i,j)...
 + awB{k}(i,j)*wStar{k-1}(i,j) - dx*dy*(pStar{k+1}(i,j)-pStar{k}(i,j)))...
 + (1-alphaU)*wOld{k}(i,j);
 end
 end
 end
 end

end

