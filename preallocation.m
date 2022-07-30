function [u,uStar,uPrime,v,vStar,vPrime,w,wStar,wPrime,p,pStar,pPrime,dU,dV,dW,uOld,vOld,wOld]...
    = preallocation(Nx,Ny,Nz)

% preallocate memory for u, v, w, p matrices
u=cell(1,Nz+2);
for m=1:Nz+2
    u{m}=zeros(Nx+1,Ny+2);
end
uStar = u;
uPrime = u;

v=cell(1,Nz+2);
for m=1:Nz+2
    v{m}=zeros(Nx+2,Ny+1);
end
vStar = v;
vPrime = v;

w=cell(1,Nz+1);
for m=1:Nz+1
    w{m}=zeros(Nx+2,Ny+2);
end
wStar = w;
wPrime = w;

p=cell(1,Nz+2);
for m=1:Nz+2
    p{m}=zeros(Nx+2,Ny+2);
end
pStar = p;
pPrime = p;

dU = u; % just for defining the dimensions of dU tensor
dV = v; % just for defining the dimensions of dV tensor
dW = w; % just for defining the dimensions of dW tensor

% initial guess ( velocity = (0 0 0) and presure = 0 )
uOld = u;
vOld = v;
wOld = w;
pStar = p;
end

