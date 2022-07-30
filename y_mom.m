function [avP,dV,vStar] = y_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,vStar,dV,Fr)

    % STEP 1b: solve y-momentum as vStar
 % setup coefficients
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny
 % convective flux using central differencing
 vCe = dy*dz*0.5*(uOld{k}(i,j)+uOld{k}(i,j+1));
 vCw = dy*dz*0.5*(uOld{k}(i-1,j)+uOld{k}(i-1,j+1));
 vCn = dx*dz*0.5*(vOld{k}(i,j)+vOld{k}(i,j+1));
 vCs = dx*dz*0.5*(vOld{k}(i,j)+vOld{k}(i,j-1));
 vCf = dx*dy*0.5*(wOld{k}(i,j)+wOld{k}(i,j+1));
 vCb = dx*dy*0.5*(wOld{k}(i-1,j)+wOld{k}(i-1,j+1));
 vDx = (dy*dz)/(dx*Re);
 vDy = (dx*dz)/(dy*Re);
 vDz = (dx*dy)/(dz*Re);
 % boundary coefficients using hybrid scheme
 avE_h = [-vCe,(vDx-0.5*vCe),0];
 avW_h = [vCw,(vDx+0.5*vCw),0];
 avN_h = [-vCn,(vDy-0.5*vCn),0];
 avS_h = [vCs,(vDy+0.5*vCs),0];
 avF_h = [-vCf,(vDz-0.5*vCf),0];
 avB_h = [vCb,(vDz+0.5*vCb),0];
 avE{k}(i,j) = max(avE_h);
 avW{k}(i,j) = max(avW_h);
 avN{k}(i,j) = max(avN_h);
 avS{k}(i,j) = max(avS_h);
 avF{k}(i,j) = max(avF_h);
 avB{k}(i,j) = max(avB_h);
 % central coefficient 
 avP{k}(i,j) = (avE{k}(i,j)+avW{k}(i,j)+avN{k}(i,j)+avS{k}(i,j)+avF{k}(i,j)+avB{k}(i,j));
 
 % set u velocity component of PCE
 dV{k}(i,j) = (dx*dz)/avP{k}(i,j);
 end
 end
 end
 % use previous calculation as initial guess in numerical method
 vStar = vOld;
 for iter = maxNMiter
 for k = 2:Nz+1
 for i = 2:Nx+1
 for j = 2:Ny
 vStar{k}(i,j) =(alphaU/avP{k}(i,j))*(avE{k}(i,j)*vStar{k}(i+1,j)...
     + avW{k}(i,j)*vStar{k}(i-1,j) + avN{k}(i,j)*vStar{k}(i,j+1)...
     +avS{k}(i,j)*vStar{k}(i,j-1) + avF{k}(i,j)*vStar{k+1}(i,j)...
     +avB{k}(i,j)*vStar{k-1}(i,j) - dx*dz*(pStar{k}(i,j+1)-pStar{k}(i,j))...
     - (dx*dy*dz)/(Fr^2))+ (1-alphaU)*vOld{k}(i,j);
 end
 end
 end
 end

end

