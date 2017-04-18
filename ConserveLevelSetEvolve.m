function  newPhi = ConserveLevelSetEvolve(phi,U, V, nx,ny,dx,dy,dt,epsilon,k,reinitialize)

newPhi = phi;
phi_s = phi;
phi_ss = phi;
lambda = dt/dx;

%% Compute the normal vector
gradPhiNorm = sqrt((phi(3:nx,2:ny-1) - phi(1:nx-2,2:ny-1)).^2 + ...
    (phi(2:nx-1,3:ny) - phi(2:nx-1,1:ny-2)).^2);%/(2dx)

Normal_x = zeros(nx,ny); Normal_y = zeros(nx,ny);
Normal_x(2:nx-1,2:ny-1) = (phi(3:nx,2:ny-1) - phi(1:nx-2,2:ny-1))./gradPhiNorm;
Normal_y(2:nx-1,2:ny-1) = (phi(2:nx-1,3:ny) - phi(2:nx-1,1:ny-2))./gradPhiNorm;




%% find the steady state to level set 

TOL =1e-3;
MaxIter = 300;
iter = 0;
difference = 1;

%% reinitialization step

if mod(k,reinitialize) == 0
while (difference > TOL && iter < MaxIter)
% for i = 1:10
    iter = iter +1;
    
    % compute the fluxes
    f = phi.*(1 - phi).* Normal_x;
    g = phi.*(1 - phi).* Normal_y;

    F_right = zeros(nx,ny); F_left = zeros(nx,ny);
    F_right(2:nx-1,2:ny-1) = (f(2:nx-1,2:ny-1) + f(3:nx,2:ny-1))/2 - ...
        epsilon * (phi(3:nx, 2:nx-1) - phi(2:nx-1,2:ny-1))/dx;
    F_left(2:nx-1,2:ny-1) = (f(2:nx-1,2:ny-1) + f(1:nx-2,2:ny-1))/2 - ...
        epsilon * (phi(2:nx-1, 2:ny-1) - phi(1:nx-2,2:ny-1))/dx;

    G_right = zeros(nx,ny); G_left = zeros(nx,ny);
    G_right(2:nx-1,2:ny-1) = (g(2:nx-1,2:ny-1) + g(2:nx-1,3:ny))/2 - ...
        epsilon * (phi(2:nx-1, 3:nx) - phi(2:nx-1,2:ny-1))/dx;
    G_left(2:nx-1,2:ny-1) = (g(2:nx-1,2:ny-1) + g(2:nx-1,1:ny-2))/2 - ...
        epsilon * (phi(2:nx-1, 2:ny-1) - phi(2:nx-1,1:ny-2))/dx;
    
    % evolve in artificial time
    newPhi = phi -lambda * (F_right - F_left + G_right - G_left);
    difference = norm(phi - newPhi);
    phi = newPhi;
end
display(sprintf('Reconstructing level set function, time = %f, difference:%g, number of iterations: %d'...
    , k*dt, difference, iter) )
end

%% Level Set evolution

for i = 2:nx-1
    for j = 2:ny-1
        
        % level set function reconstruction using Superbee limiter
        if i == 2
            Sx_left = Lim((phi(i,j)-phi(i-1,j))/dx,(phi(i-1,j)-0)/dx);
        else
            Sx_left = Lim((phi(i,j)-phi(i-1,j))/dx,(phi(i-1,j)-phi(i-2,j))/dx);
        end
        if i == nx-1
            Sx_right = Lim((0-phi(i+1,j))/dx,(phi(i+1,j)-phi(i,j))/dx);
        else
            Sx_right = Lim((phi(i+2,j)-phi(i+1,j))/dx,(phi(i+1,j)-phi(i,j))/dx);
        end
        Sx = Lim((phi(i+1,j)-phi(i,j))/dx,(phi(i,j)-phi(i-1,j))/dx);
        
        if j == 2
            Sy_bot = Lim((phi(i,j)-phi(i,j-1))/dy,(phi(i,j-1)-0)/dy);
        else
            Sy_bot = Lim((phi(i,j)-phi(i,j-1))/dy,(phi(i,j-1)-phi(i,j-2))/dy);
        end
        if j == nx-1
            Sy_top = Lim((0-phi(i,j+1))/dy,(phi(i,j+1)-phi(i,j))/dy);
        else
            Sy_top = Lim((phi(i,j+2)-phi(i,j+1))/dy,(phi(i,j+1)-phi(i,j))/dy);
        end
        Sy = Lim((phi(i,j+1)-phi(i,j))/dy,(phi(i,j)-phi(i,j-1))/dy);

        % reconstruct phi
        phi_left_plus = phi(i,j) - Sx * dx/2;
        phi_left_minus = phi(i-1,j) + Sx_left * dx/2;
        
        phi_right_plus = phi(i+1,j) - Sx_right * dx/2;
        phi_right_minus = phi(i,j) + Sx * dx/2; 
        
        phi_top_plus = phi(i,j+1) - Sy_top * dy/2;
        phi_top_minus = phi(i,j) + Sy * dy/2;
        
        phi_bot_plus = phi(i,j) - Sy * dy/2;
        phi_bot_minus = phi(i,j-1) + Sy_bot * dy/2;
        
        % compute the fluxes
        F_right = max(U(i,j),0) * phi_right_minus + ...
            min(U(i,j),0) * phi_right_plus;
        F_left = max(U(i-1,j),0) * phi_left_minus + ...
            min(U(i-1,j),0) * phi_left_plus;
        
        G_top = max(V(i,j),0) * phi_top_minus + ...
            min(V(i,j),0) * phi_top_plus;
        G_bot = max(V(i,j-1),0) * phi_bot_minus + ...
            min(V(i,j-1),0) *phi_bot_plus;
        
        % compute the total flux
        F = -(F_right-F_left + G_top - G_bot)/dx;
        
        % perform Runge-Kuttah integration
        phi_s(i,j) = phi(i,j) + dt * F;
        
    end
end


for i = 2:nx-1
    for j = 2:ny-1
        
        % level set function reconstruction using Superbee limiter
        if i == 2
            Sx_left = Lim((phi_s(i,j)-phi_s(i-1,j))/dx,(phi_s(i-1,j)-0)/dx);
        else
        Sx_left = Lim((phi_s(i,j)-phi_s(i-1,j))/dx,(phi_s(i-1,j)-phi_s(i-2,j))/dx);
        end
        if i == nx-1
            Sx_right = Lim((0-phi_s(i+1,j))/dx,(phi_s(i+1,j)-phi_s(i,j))/dx);
        else
            Sx_right = Lim((phi_s(i+2,j)-phi_s(i+1,j))/dx,(phi_s(i+1,j)-phi_s(i,j))/dx);
        end
        Sx = Lim((phi_s(i+1,j)-phi_s(i,j))/dx,(phi_s(i,j)-phi_s(i-1,j))/dx);
        
        if j == 2
            Sy_bot = Lim((phi_s(i,j)-phi_s(i,j-1))/dy,(phi_s(i,j-1)-0)/dy);
        else
            Sy_bot = Lim((phi_s(i,j)-phi_s(i,j-1))/dy,(phi_s(i,j-1)-phi_s(i,j-2))/dy);
        end
        if j == ny-1
            Sy_top = Lim((0-phi_s(i,j+1))/dy,(phi_s(i,j+1)-phi_s(i,j))/dy);
        else
            Sy_top = Lim((phi_s(i,j+2)-phi_s(i,j+1))/dy,(phi_s(i,j+1)-phi_s(i,j))/dy);
        end
        Sy = Lim((phi_s(i,j+1)-phi_s(i,j))/dy,(phi_s(i,j)-phi_s(i,j-1))/dy);

        % reconstruct phi_s
        phi_s_left_plus = phi_s(i,j) - Sx * dx/2;
        phi_s_left_minus = phi_s(i-1,j) + Sx_left * dx/2;
        
        phi_s_right_plus = phi_s(i+1,j) - Sx_right * dx/2;
        phi_s_right_minus = phi_s(i,j) + Sx * dx/2; 
        
        phi_s_top_plus = phi_s(i,j+1) - Sy_top * dy/2;
        phi_s_top_minus = phi_s(i,j) + Sy * dy/2;
        
        phi_s_bot_plus = phi_s(i,j) - Sy * dy/2;
        phi_s_bot_minus = phi_s(i,j-1) + Sy_bot * dy/2;
        
        % compute the fluxes
        F_right = max(U(i,j),0) * phi_s_right_minus + ...
            min(U(i,j),0) * phi_s_right_plus;
        F_left = max(U(i-1,j),0) * phi_s_left_minus + ...
            min(U(i-1,j),0) * phi_s_left_plus;
        
        G_top = max(V(i,j),0) * phi_s_top_minus + ...
            min(V(i,j),0) * phi_s_top_plus;
        G_bot = max(V(i,j-1),0) * phi_s_bot_minus + ...
            min(V(i,j-1),0) *phi_s_bot_plus;
        
        % compute the total flux
        F = -(F_right-F_left + G_top - G_bot)/dx;
        
        % perform Runge-Kuttah integration
        phi_ss(i,j) = phi_s(i,j) + dt * F;
        
    end
end

        
newPhi = 0.5 * (phi + phi_ss);


