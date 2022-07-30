function [] = postProcessing(n,telapsed,max_residual,x,y,z,Nx,Ny,Nz,...
    u,v,w,p,dx,dy,dz,L1,L2,L3)

% number of iterations and computational time
 Num_Iteration = n;
 fprintf('\n number of iterations = %05e \n',Num_Iteration);
 disp(' ************************************************ ')
 fprintf('\n max of residual = %05e \n',max_residual);
 disp(' ************************************************ ')
 fprintf('\n computational time = %05e',telapsed/60);
 disp(' minutes ')
 disp(' ************************************************ ')
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PostProcessing

% streamline plot
[xM,yM] = meshgrid(x,y);
U = u{floor((Nz+1)/2)}(1:Nx+1,2:Ny+2);
V = v{floor((Nz+1)/2)}(1:Nx+1,1:Ny+1);
W = w{floor((Nz+1)/2)}(1:Nx+1,1:Ny+1);
P = p{floor((Nz+1)/2)}(1:Nx+1,1:Ny+1);
sx = 0:dx:L1;
sy = 0:dy:L2;
figure(1)
axis equal;
grid on
xlim([0 L1]);
ylim([0 L2]);
h = streamline(xM,yM,U',V',sx,sy);
set(h,'Color','b');
xlabel('x [m]'); ylabel('y [m]');
title( ' streamlines at x-y plane at z = 0.5 m ')

% computational Mesh
figure(2)
surf(xM,yM,ones(Nx+1,Ny+1));
view(2); axis equal; axis ([0 L1 0 L2]);
xlabel('x [m]'); ylabel('y [m]');
title( ' computatiomnal grid at x-y plane at z = 0.5 m ')

% contour plot for velocity magnitude
Umag=sqrt(U'.*U' + V'.*V' + W'.*W');
figure(3)
axis equal;
grid on
contourf(xM,yM,Umag,'--','ShowText','on'),colormap(jet);
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' Velocity contour at x-y plane at z = 0.5 m ')

% contour plot for velocity magnitude
figure(4)
pcolor(xM,yM,Umag),colormap(jet),shading interp;
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' Velocity contour at x-y plane at z = 0.5 m ')

% contour plot for u velocity at z=0.5m
figure(5)
pcolor(xM,yM,U'),colormap(jet),shading interp;
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' V_x = u contour at x-y plane at z = 0.5m ')

% contour plot for v velocity
figure(6)
pcolor(xM,yM,V'),colormap(jet),shading interp;
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' V_y = v contour at x-y plane at z = 0.5m ')

% contour plot for w velocity
figure(7)
pcolor(xM,yM,W'),colormap(jet),shading interp;
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' V_z = w contour at x-y plane at z = 0.5m ')

% contour plot for pressure
figure(8)
pcolor(xM,yM,P'),colormap(jet),shading interp;
colorbar,pbaspect([1/L1,1/L2,1]);
xlabel('x [m]'); ylabel('y [m]');
title( ' pressure contour at x-y plane at z = 0.5m ')


%velocity profile at x=(L/2)
figure(9)   
analytical_solution=csvread('Re400.csv');
ya=analytical_solution(:,1);
ua=analytical_solution(:,2);
plot(ya,ua,'*','DisplayName','theoretical')
xlabel('u (m/s)')
ylabel('y (m) ')
title('velocity profile at x=z=0.5 m')
hold on
yyy=linspace(0,L2,Nx+2);
uuu=u{floor((Nz+1)/2)}(floor((Nx+1)/2),:);
plot(uuu,yyy,'DisplayName','numerical')
legend('analytical solution','numerical solution');
grid on

end

