function [ deformed_vol ] = deform(vol, Vmax)
    [Nx, Ny, Nz] = size(vol);
    [X, Y, Z] = meshgrid(1:Nx, 1:Ny, 1:Nz);
    fx = 2*pi/(Nx-1);
    fy = 2*pi/(Ny-1);
    fz = 2*pi/(Nz-1);
    A = rand(3,1);
    V = Vmax*rand(1);
    phi = 2*pi*rand(3,1);
    Vx = V*(A(1)/norm(A))*sin(fx*(X-1)+phi(1));
    Vy = V*(A(2)/norm(A))*sin(fy*(Y-1)+phi(2));
    Vz = V*(A(3)/norm(A))*sin(fz*(Z-1)+phi(3));
    deformed_vol = interp3(vol,X+Vy,Y+Vz,Z+Vx,'spline',0);
end

