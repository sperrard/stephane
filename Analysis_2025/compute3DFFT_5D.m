function [kx, ky, kz, C_fft] = compute3DFFT_5D(X, Y, Z, C)
    % Compute the 3D Fourier Transform of C with proper frequency axis calculation
    % and fftshift application.
    %
    % Inputs:
    %   X, Y, Z - 1D arrays representing the spatial coordinates.
    %   C - 3D array of size (length(X), length(Y), length(Z))
    %
    % Outputs:
    %   Fx, Fy, Fz - Frequency axes corresponding to X, Y, and Z.
    %   C_fft - Shifted 3D Fourier transform of C.

    % Ensure inputs are column vectors
    X = X(:);
    Y = Y(:);
    Z = Z(:);
    
    % Get number of points along each axis
    Nx = length(X);
    Ny = length(Y);
    Nz = length(Z);
    
    % Compute spatial resolutions (assuming uniform spacing)
    dx = mean(diff(X));
    dy = mean(diff(Y));
    dz = mean(diff(Z));
    
    % Compute frequency axes
    kx = 2*pi*(-floor(Nx/2):ceil(Nx/2)-1) / (Nx*dx);
    ky = 2*pi*(-floor(Ny/2):ceil(Ny/2)-1) / (Ny*dy);
    kz = 2*pi*(-floor(Nz/2):ceil(Nz/2)-1) / (Nz*dz);
    
    % Compute the 3D FFT and apply fftshift
    [nx,ny,nz,nt,NFilm] = size(C);
    C_fft = zeros(nx,ny,nz,nt,NFilm);
    for i=1:nt
        for j=1:NFilm
            Cp = C(:,:,:,i,j);
            Cp_fft = fftshift(fftn(Cp));
            C_fft(:,:,:,i,j) = Cp_fft/numel(Cp_fft);
        end
    end

end
