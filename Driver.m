% 
% Driver.m
% Initializes the level set problem on  a rectangular grid

%%
%Grid construction
lx = 4;
ly = 4;
nx = 100; ny = 100;
dx = lx/nx; dy = ly/ny;
dt = 0.01;
x = -lx/2:dx:lx/2-dx;
y = -ly/2:dy:ly/2-dy;
[X,Y] = meshgrid(x,y);

phi = sqrt((X-1).^2+(Y-1).^2)-0.2;

contour(X,Y,phi,'showtext','on')
axis equal

% u = X;
% v = Y;
u = -Y;
v = X;
Tfinal = 2*pi;
quiver(X,Y,u,v)
%G = numgrid('S',nx)
%newphi
newPhi = LevelSetEvolve(phi,v,u,nx,ny,dx,dy,dt);
contour(X,Y,phi,[0,0],'linewidth',2)
xlabel('X')
ylabel('Y')
axis equal
N = ceil(Tfinal/dt);
for k=1:N
    time = k*dt;
    newPhi = LevelSetEvolve(phi,v,u,nx,ny,dx,dy,dt);
    phi = newPhi;
%     contourf(X,Y,phi)
    colorbar
    contour(X,Y,phi,[0,0])
    xlabel('X')
    ylabel('Y')
%     contour(X,Y,phi,'showtext','on')
    axis equal
    title(sprintf('time=%d',time))
    pause(0.001)
end
