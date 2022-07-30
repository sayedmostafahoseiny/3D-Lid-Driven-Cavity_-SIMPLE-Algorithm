function [L1,L2,L3,ulid,g,Re,Fr,alphaU,alphaP,maxNMiter,err_criteria,Min_Iteration,max_residual]...
    = parameters(nu)

% variable declaration
L = 1;  
L1 = L; % Length in x direction (m)
L2 = L; % Length in y direction ( Height ) (m)
L3 = L; % Length in z direction (m)
ulid = 1; % lid velocity or upper wall velocity in (m/s)
g = 9.81; % gravitational acceleration in (m/s^2) ( in y direction )
Re = (ulid * L)/nu; % Reynolds number
Fr = (ulid)/sqrt(g*L); % Froude number

alphaU = 0.5; % Velocity under relaxation factor
alphaP = 0.8; % Pressure under relaxation factor
maxNMiter = 1; % Iteration for numerical method
err_criteria = 1e-5; % Convergence criteria
Min_Iteration = 50; % the minimum Iteration of SIMPLE loop
% It means SIMPLE loop iterates at least Min_Iteration times
max_residual = 1 ; % it's just an initial value (It will be calculated in any iteration)
% (it must be greater than the err_criteria to start the loop)

end

