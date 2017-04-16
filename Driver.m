% 
% Driver.m
% Initializes the level set problem on  a rectangular grid

%%
%Grid construction
clear all

lx = 2;
ly = 2;
nx = 100; ny = 100;
dx = lx/nx; dy = ly/ny;
dt = 0.01;
epsilon = dx/2;
x = -lx/2:dx:lx/2-dx;
y = -ly/2:dy:ly/2-dy;
[X,Y] = meshgrid(x,y);

phi = sqrt((X-0.5).^2+(Y-0.5).^2)-0.3;
phi0 = 1./(1+exp(phi./epsilon));
phi = phi0;
% phi = phi_sd>0;
contour(X,Y,phi,[0.5 0.5],'showtext','on')
axis equal
hold on

xflow = -lx/2 + dx/2:dx:lx/2 - 3*dx/2;
yflow = -ly/2 + dy/2:dy:ly/2 - 3*dy/2;
[Xflow,Yflow] = meshgrid(xflow,yflow);
% u = X;
% v = Y;
U = -Yflow;
V = Xflow;

U = sin(pi*Xflow).^2.*sin(2*pi*Yflow);
V = -sin(pi*Yflow).^2.*sin(2*pi*Xflow);
% Tfinal = 2*pi;
Tfinal = 1;
quiver(Xflow,Yflow,U,V)
%G = numgrid('S',nx)
% %newphi
% newPhi = LevelSetEvolve(phi,v,u,nx,ny,dx,dy,dt);
% contour(X,Y,phi,[0,0],'linewidth',2)
% xlabel('X')
% ylabel('Y')
% axis equal
% title('time=0')
% axis([-2 2 -2 2])
N = ceil(Tfinal/dt);
reinitialize = 5;
hold off
N=100;
area = zeros(1,N);

for k=1:N
    time = k*dt;
%     newPhi = LevelSetEvolve(phi,v,u,nx,ny,dx,dy,dt);
    newPhi = ConserveLevelSetEvolve(phi,V,U,nx,ny,dx,dy,dt,epsilon,k,reinitialize);
    phi = newPhi;
%     contourf(X,Y,phi)
    colorbar
    [c,h] = contour(X,Y,phi,[0.5 0.5],'showtext','on');
    xlabel('X')
    ylabel('Y')
%     contour(X,Y,phi,'showtext','on')
    axis equal
    axis([0 1 0 1])
    title(sprintf('time=%d',time))
    drawnow
    area(k) = polyarea(c(1,:),c(2,:));
end
hold on
contour(X,Y,phi0,[0.5 0.5],'showtext','on')
axis equal
quiver(Xflow,Yflow,U,V)
%% run conservative evolution

lx = 8;
ly = 8;
nx = 100; ny = 100;
dx = lx/nx; dy = ly/ny;

x = -lx/2:dx:lx/2;
y = -ly/2:dy:ly/2;
[X,Y] = meshgrid(x,y);

dt = 0.0001;
epsilon = 0.1;

%% initialize level set function

phi_sd = sqrt((X-1).^2+(Y-1).^2)-0.3;

