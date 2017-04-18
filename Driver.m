% 
% Driver.m
% Initializes the level set problem on  a rectangular grid

%%
%Grid construction
clear all

lx = 4;
ly = 4;
nx = 100; ny = 100;
dx = lx/nx; dy = ly/ny;
dt = 0.01;
epsilon = dx/2;
x = -lx/2:dx:lx/2-dx;
y = -ly/2:dy:ly/2-dy;
[X,Y] = meshgrid(x,y);

phi = sqrt((X-0.5).^2+(Y-0.75).^2)-0.15;
phi0 = 1./(1+exp(phi./epsilon));
phi = phi0;

contour(X,Y,phi,[0.5 0.5],'showtext','on')
axis equal
hold on

xflow = -lx/2 + dx/2:dx:lx/2 - 3*dx/2;
yflow = -ly/2 + dy/2:dy:ly/2 - 3*dy/2;
[Xflow,Yflow] = meshgrid(xflow,yflow);

% Rotating circle
U = -Yflow;
V = Xflow;
Tfinal = 2*pi;
% Tfinal = 3*Tfinal/4;

% Vortex flow
% U = sin(pi*Xflow).^2.*sin(2*pi*Yflow);
% V = -sin(pi*Yflow).^2.*sin(2*pi*Xflow);
% Tfinal = 1;
% Tfinal = 3*Tfinal/4;
quiver(Xflow,Yflow,U,V)
N = ceil(Tfinal/dt);
reinitialize = 5;
hold off
%N=100;
area = zeros(1,N);

%% Conservative level set
for k=1:N
    time = k*dt;
%     newPhi = LevelSetEvolve(phi,V,U,nx,ny,dx,dy,dt);
    newPhi = ConserveLevelSetEvolve(phi,V,U,nx,ny,dx,dy,dt,epsilon,k,reinitialize);
    phi = newPhi;
    
    colorbar
    [c,h] = contour(X,Y,phi,[0.5 0.5],'showtext','on','linewidth',3);
    xlabel('X')
    ylabel('Y')
    axis equal
    axis([-lx/2 lx/2 -ly/2 ly/2])
    title(sprintf('time=%d',time))
    hold on
    quiver(Xflow,Yflow,U,V)
    set(gca,'fontsize',14)
    hold off
    drawnow
    area(k) = polyarea(c(1,:),c(2,:));
end
hold on
contour(X,Y,phi0,[0.5 0.5],'showtext','on','linewidth',3)
axis equal
quiver(Xflow,Yflow,U,V)
axis([-lx/2 lx/2 -ly/2 ly/2])

% plot area 
figure(2)
title('Volume conservation vortex', 'fontsize',13)
hold on
plot(dt:dt:N*dt,area,'linewidth',3)
xlabel('time','fontsize',13)
ylabel('area','fontsize',13)
axis([0 Tfinal 0 1.1*max(area)])
set(gca,'fontsize',13)