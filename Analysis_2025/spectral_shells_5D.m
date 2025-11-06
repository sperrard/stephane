function [Eshell, Emod] = spectral_shells_5D(kx, ky, kz, E, N, kmax, space_type)
% Compute the spectral shell decomposition of a 3D energy spectrum
%
% Inputs:
%   kx, ky, kz - 1D arrays representing the wave number components
%   E - 3D array of size (length(kx), length(ky), length(kz))
%   N - Number of bins for kmod_1D
%   klim (optional) - 2-element array specifying the range of kmod_1D
%   space_type (optional) - 'lin' for linear spacing, 'log' for logarithmic spacing
%
% Outputs:
%   Eshell - Cell array containing spectral shell data
%   Emod - Struct with spectral shell-averaged values

% Ensure kx, ky, kz are column vectors
kx = kx(:);
ky = ky(:);
kz = kz(:);

% Create 3D meshgrids for kx, ky, kz
[KX, KY, KZ] = ndgrid(kx, ky, kz);

% Compute the modulus of the wave vectors
kmod_3D = sqrt(KX.^2 + KY.^2 + KZ.^2);

% Determine kmod_1D range
k_nonzero = kmod_3D(kmod_3D > 0);
kmod_min = min(k_nonzero);
if nargin < 6 || isempty(kmax)
    kmod_max = max(kmod_3D(:));
else
    kmod_max = kmax;
end

% Default space_type to 'lin' if not provided
if nargin < 7 || isempty(space_type)
    space_type = 'lin';
end

% Define kmod_1D based on spacing type
if strcmp(space_type, 'log')
    kmod_1D = logspace(log10(kmod_min), log10(kmod_max), N);
else
    kmod_1D = linspace(kmod_min, kmod_max, N);
end

% Compute the mean difference dkmod
dkmod = gradient(kmod_1D);

% Initialize cell array Eshell
Eshell = {};

% Initialize structure Emod
Emod = struct('kmod', [], 'spec', []);

c = 1;
% Iterate over kmod_1D values
[nx,ny,nz,nt,NFilm] = size(E);

nc = length(kmod_1D);

Emod_mat = zeros(nc,nt,NFilm);
for j=1:nt
    disp(j)
    for k=1:NFilm
        Ep = E(:,:,:,j,k);
        for i = 1:nc
            % Find indices where kmod_3D is within the current shell range
            mask = (kmod_3D >= kmod_1D(i) - dkmod(i)/2) & (kmod_3D <= kmod_1D(i) + dkmod(i)/2);

            if sum(sum(sum(mask))) == 0
                continue;
            end

            % Store the corresponding values in Eshell
            Eshell{c}.kx = KX(mask);
            Eshell{c}.ky = KY(mask);
            Eshell{c}.kz = KZ(mask);
            Eshell{c}.kmod = kmod_3D(mask);
            Eshell{c}.spec = Ep(mask);

            % Compute the mean spectral energy in the shell for Emod
            Emod.kmod(c) = kmod_1D(i);%kind of useless, but who cares ?
            Emod_mat(c,j,k) = 4*pi*(kmod_1D(i)^2)*abs(mean(Eshell{c}.spec));

            c = c + 1;
        end
    end
end
Emod.mat=nanmean(Emod_mat,2);
end
