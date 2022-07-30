% solving the Steady, Incompressible, Isothermal, 3D, Laminar flow in Lid driven cavity
% The gravitational acceleration is in y dirrection
% finite volume method + SIMPLE Algorithm with Hybrid Scheme For Convective and diffusive fluxes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc; format long

% Defining parameters
nu = input ( ' Enter the kinematic viscisity (m^2/s) = ' ) ;
[L1,L2,L3,ulid,g,Re,Fr,alphaU,alphaP,maxNMiter,err_criteria,Min_Iteration,max_residual]...
    = parameters(nu);
fprintf('\n Reynolds number = %05e \n',Re); 
disp ( ' ********************************************* ')

% Defining computational Grid (Mesh)
[Nx,Ny,Nz,dx,x,dy,y,dz,z] = Mesh(L1,L2,L3);

% preallocation of u,v,w,p and their star and their prime Tensors
[u,uStar,uPrime,v,vStar,vPrime,w,wStar,wPrime,p,pStar,pPrime,dU,dV,dW,uOld,vOld,wOld]...
    = preallocation(Nx,Ny,Nz);

% setting Boundary conditions
[u,v,w] = setting_BCs(u,v,w,Nx,Ny,Nz,ulid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tstart = tic; % for calculation of computational time
% SIMPLE algorithm ( SIMPLE Loop )
n = 1 ; % counter for the SIMPLE loop
while ( max_residual > err_criteria || n < Min_Iteration)

    % STEP 1a: solve x-momentum as uStar ( prediction step )
[auP,dU,uStar] = x_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,uStar,dU);
    
    % STEP 1b: solve y-momentum as vStar ( prediction step )
[avP,dV,vStar] = y_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,vStar,dV,Fr);

    % STEP 1c: solve z-momentum as wStar ( prediction step )
[awP,dW,wStar] = z_mom(Nx,Ny,Nz,dx,dy,dz,uOld,vOld,wOld,Re,alphaU,maxNMiter,pStar,wStar,dW);
 % ************************************************************
 
 % STEP 2: solve pressure correction equation (PCE = pressure correction equation)
 [pPrime] = PCE(Nx,Ny,Nz,dx,dy,dz,dU,dV,dW,uStar,vStar,wStar,pPrime,n,maxNMiter);

 % STEP 3: calculate velocity corrections by pPrime
 [uPrime,vPrime,wPrime] = Velocity_correctors(Nx,Ny,Nz,uPrime,dU,pPrime,vPrime,dV,wPrime,dW);

 % STEP 4: correct u, v, w and p ( correction step )
 [u,v,w,p] = var_corrections(Nx,Ny,Nz,alphaP,alphaU,p,pStar,...
    pPrime,u,uStar,uPrime,v,vStar,vPrime,w,wStar,wPrime);
 
 % check for convergence ( max residual for breaking the SIMPLE loop )
[max_residual] = max_residual_calculation(Nx,Ny,Nz,u,uOld,v,vOld,w,wOld);

 % STEP 5: update cell velocities ( update old values )
 uOld = u;
 vOld = v;
 wOld = w;
 pStar = p;
 n = n + 1 ; % update the counter
 
end % End of the SIMPLE loop
telapsed = toc(tstart); % for the calculation of computational time

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% post processing ( profiles - contours - streamlines )
postProcessing(n,telapsed,max_residual,x,y,z,Nx,Ny,Nz,...
    u,v,w,p,dx,dy,dz,L1,L2,L3);






