function [u,v,w] = setting_BCs(u,v,w,Nx,Ny,Nz,ulid)

 % Lid Driven Cavity Boundary conditions
 % apply boundary conditions
 for k = 1:Nz+2
 u{k}(:,1) = 0.0; % bottom boundary
 v{k}(:,1) = 0.0;
 w{k}(:,1) = 0.0;
 end
 for k = 1:Nz+2
 u{k}(:,Ny+2) = ulid; % top boundary
 v{k}(:,Ny+1) = 0.0;
 w{k}(:,Ny+2) = 0.0;
 end
 for k = 1:Nz+2
 u{k}(1,:) = 0.0; % left boundary
 v{k}(1,:) = 0.0;
 w{k}(1,:) = 0.0;
 end
 for k = 1:Nz+2
 u{k}(Nx+1,:) = 0.0; % right boundary
 v{k}(Nx+2,:) = 0.0;
 w{k}(Nx+2,:) = 0.0;
 end
 % front wall boundary
 u{1}(:,:) = 0.0; 
 v{1}(:,:) = 0.0;
 w{1}(:,:) = 0.0;
 % back wall boundary
 u{Nz+2}(:,:) = 0.0; 
 v{Nz+2}(:,:) = 0.0;
 w{Nz+1}(:,:) = 0.0;
 
end

