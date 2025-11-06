N = 100;
type = 'log';

kmax_hp = 1;
kmax_vp = 1;

geometry = 'mpmp';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% Horizontal Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('C3D_SP_2021_0623.mat');


% Integral length scale - x
Gamma_xx_hp = Mhp.Cxx(:,39,1);
dx_hp = Mhp.dx;

xft_xx = linspace(0, 50, 100);
yft_xx = exp(-0.0217)*exp(-xft_xx/8.62);

figure(); 
plot(dx_hp, Gamma_xx_hp, '--rs');
hold on
plot(-dx_hp, Gamma_xx_hp, '-k+');

hold on;
plot(xft_xx, yft_xx, '--k', LineWidth=2);
% xlim([-50, 50])
legend('C_{xx}', '\propto exp(-x/8.62)')
yscale log;
title('hp')



% Integral length scale - y
Gamma_yy_hp = Mhp.Cyy(37,:,1);
dy_hp = Mhp.dy;

xft_yy = linspace(0, 45, 100);
yft_yy = exp(0.0319)*exp(-xft_xx/8.57);% do a proper fit !

figure(); 
plot(dy_hp, Gamma_yy_hp, '--rs');
hold on;
plot(-dy_hp, Gamma_yy_hp, '--k+');

plot(xft_xx, yft_xx, '--k', LineWidth=2);
% xlim([-45, 45])
legend('C_{yy}', '\propto exp(-x/8.57)')
yscale log;
title('hp')

%%

% Make symmetric along z
Mhp2 = Mhp;
Mhp2.Cxx = cat(3, flip(Mhp.Cxx, 3), Mhp.Cxx);
Mhp2.Cyy = cat(3, flip(Mhp.Cyy, 3), Mhp.Cyy);
Mhp2.Cxy = cat(3, flip(Mhp.Cxy, 3), Mhp.Cxy);
Mhp2.Cyx = cat(3, flip(Mhp.Cyx, 3), Mhp.Cyx);
Mhp2.dz = cat(2, -flip(Mhp.dz, 2), Mhp.dz);


Mhp2.Cxx_sym = (Mhp2.Cxx+flip(flip(flip(Mhp2.Cxx,1),2),3))/2;
Mhp2.Cyy_sym = (Mhp2.Cyy+flip(flip(flip(Mhp2.Cyy,1),2),3))/2;


%we impose : 
Mhp2.Cxy_sym = (Mhp2.Cxy-flip(flip(flip(Mhp2.Cyx,1),2),3))/2;%+Czz(-x,-y,-z))/2;
Mhp2.Cyx_sym = (Mhp2.Cyx-flip(flip(flip(Mhp2.Cxy,1),2),3))/2;%+Czz(-x,-y,-z))/2;

% 
% figure()
% plot(dx_hp,Mhp2.Cxx_sym(:,39,126))
% hold on
% plot(-dx_hp,Mhp2.Cxx_sym(:,39,126))
% yscale log
% 
% figure()
% plot(dy_hp,abs(Mhp2.Cxy_sym(37,:,126)))
% hold on
% plot(dx_hp,abs(Mhp2.Cxy(:,39,126)))
% yscale log

%%
% Energy spectrum : xx
disp('Compute Cxx')
[kx_xx_hp, ky_xx_hp, kz_xx_hp, E_xx_hp] = compute3DFFT_5D(Mhp.dx, Mhp.dy, Mhp.dz, Mhp.Cxx);

% Energy spectrum : yy
disp('Compute Cyy')
[kx_hp, ky_hp, kz_hp, E_yy_hp] = compute3DFFT_5D(Mhp.dx, Mhp.dy, Mhp.dz, Mhp.Cyy);

% Energy spectrum : xy
disp('Compute Cxy')
[~, ~, ~, E_xy_hp] = compute3DFFT_5D(Mhp.dx, Mhp.dy, Mhp.dz, Mhp.Cxy);

% Energy spectrum : yx
disp('Compute Cyx')
[~, ~, ~, E_yx_hp] = compute3DFFT_5D(Mhp.dx, Mhp.dy, Mhp.dz, Mhp.Cyx);

%%
% average over shells
[Eshell_xx_hp, Emod_xx_hp] = spectral_shells(kx_xx_hp, ky_xx_hp, kz_xx_hp, mean(mean(abs(E_xx_hp),4),5), N, kmax_hp, type);
[Eshell_yy_hp, Emod_yy_hp] = spectral_shells(kx_hp, ky_hp, kz_hp, mean(mean(E_yy_hp,4),5), N, kmax_hp, type);
[Eshell_xy_hp, Emod_xy_hp] = spectral_shells(kx_hp, ky_hp, kz_hp, mean(mean(E_xy_hp,4),5), N, kmax_hp, type);
[Eshell_yx_hp, Emod_yx_hp] = spectral_shells(kx_hp, ky_hp, kz_hp, mean(mean(E_yx_hp,4),5), N, kmax_hp, type);

%%

NFilm = 11;
figure(); hold on; 

for j=1:NFilm

plot(Emod_xx_hp.kmod, Emod_xx_hp.spec, '--rs'); 
plot(Emod_yy_hp.kmod, Emod_yy_hp.spec, '--bo');
xline(1/10, '--k', LineWidth=2)
% xline(2*pi/10, '--k', LineWidth=2)
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$E(k)$', 'Interpreter', 'latex');
legend({'$E_{xx}$', '$E_{yy}$'},'Interpreter','latex')
title('hp')

xscale log; yscale log;
end
%%
% Make xy-yx symmetric
% E_xy_hp = 0.5*(E_xy_hp + conj(E_yx_hp));
% E_yx_hp = conj(E_xy_hp);

% Evaluate remaining elements of Energy spectrum matrix
S = size(E_xx_hp);
E_xz_hp = zeros(S); E_zx_hp = zeros(S); E_yz_hp = zeros(S); E_zy_hp = zeros(S); E_zz_hp = zeros(S);
lambda_hp = cell(S);

for j = 1:S(1)
    for k = 1:S(2)
        for l = 1:S(3)
            E_xz_hp(j,k,l) = - (kx_hp(j)/kz_hp(l))*E_xx_hp(j,k,l) - (ky_hp(k)/kz_hp(l))*E_xy_hp(j,k,l) ;
            E_zx_hp(j,k,l) = conj(E_xz_hp(j,k,l));

            E_yz_hp(j,k,l) = - (kx_hp(j)/kz_hp(l))*E_yx_hp(j,k,l) - (ky_hp(k)/kz_hp(l))*E_yy_hp(j,k,l) ;
            E_zy_hp(j,k,l) = conj(E_yz_hp(j,k,l));

            E_zz_hp(j,k,l) = - (kx_hp(j)/kz_hp(l))*E_zx_hp(j,k,l) - (ky_hp(k)/kz_hp(l))*E_zy_hp(j,k,l) ;

            E_matrix = [E_xx_hp(j,k,l) E_xy_hp(j,k,l) E_xz_hp(j,k,l) ; E_yx_hp(j,k,l) E_yy_hp(j,k,l) E_yz_hp(j,k,l) ; E_zx_hp(j,k,l) E_zy_hp(j,k,l) E_zz_hp(j,k,l)];
            %E_matrix = [E_xx E_xy E_xz ; 
            %            E_yx E_yy E_yz ; 
            %            E_zx E_zy E_zz];
            eigval = eig(E_matrix);
            lambda_hp{j,k,l} = eigval;
        end
    end
end

[lambda_shell_hp, lambda_mod_hp] = analyze_lambda_ratio(kx_hp, ky_hp, kz_hp, lambda_hp, N, kmax_hp, type);


figure(); hold on; 
plot(Emod_xx_hp.kmod, Emod_xx_hp.spec, '--rs'); 
plot(Emod_yy_hp.kmod, Emod_yy_hp.spec, '--bo');
xline(1/10, '--k', LineWidth=2)
% xline(2*pi/10, '--k', LineWidth=2)
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$E(k)$', 'Interpreter', 'latex');
legend({'$E_{xx}$', '$E_{yy}$'},'Interpreter','latex')
title('hp')

xscale log; yscale log;

style = {'--kp','--mx'};
figure(); 
subplot(2,1,1); hold on;
plot(lambda_mod_hp.kmod, lambda_mod_hp.lambdaratio, style{1});
legend('\langle \lambda_{3}/\lambda_{2} \rangle')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$\langle \lambda_{3}/\lambda_{2} \rangle$', 'Interpreter', 'latex');
% yscale log;
subplot(2,1,2); hold on;
plot(lambda_mod_hp.kmod, lambda_mod_hp.lambdamin, style{2});
legend('\langle \lambda_{2} \rangle')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$\langle \lambda_{2} \rangle$', 'Interpreter', 'latex');
% yscale log;
title('hp')

getframe();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Vertical Plane %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('C3D_SP_2021_0928.mat');


% Geometry
if strcmp(geometry, 'mpmp')
    indices = 1:3;
elseif strcmp(geometry, 'pppp')
    indices = 4:5;
elseif strcmp(geometry, 'pmpm')
    indices = 6:7;
end


% Integral length scale - x
Gamma_xx_vp = [];
dx_vp = [];
for i = indices
    Gamma_xx_vp = [Gamma_xx_vp ; reshape(Mvp.Cxx(:,61,1,i), 1, [])];
    dx_vp = [dx_vp ; reshape(Mvp.dx, 1, [])];
end

% Integral length scale - z
Gamma_zz_vp = [];
dz_vp = [];
for i = indices
    Gamma_zz_vp = [Gamma_zz_vp ; reshape(Mvp.Czz(58,:,1,i), 1, [])];
    dz_vp = [dz_vp ; reshape(Mvp.dz, 1, [])];
end


xft_xx = linspace(0, 45, 100);
yft_xx = 0.647*exp(-xft_xx/10.4);
xft_zz = linspace(0, 45, 100);
yft_zz = 1.32*exp(-xft_zz/7.63);

style = {'--rs','--bo','--gd'};
figure(); 
subplot(2,1,1); hold on;
for i = 1:length(indices)
    plot(dx_vp(i,:), Gamma_xx_vp(i,:), style{i});
end
plot(xft_xx, yft_xx, '--k', LineWidth=2);
xlim([-45, 45])
legend('C_{xx}', '\propto exp(-x/10.4)')
yscale log;
subplot(2,1,2); hold on;
for i = 1:length(indices)
    plot(dz_vp(i,:), Gamma_zz_vp(i,:), style{i});
end
plot(xft_zz, yft_zz, '--k', LineWidth=2);
xlim([-45, 45])
legend('C_{zz}', '\propto exp(-x/7.64)')
yscale log;



[Eshell_xx_vp, Emod_xx_vp, Eshell_zz_vp, Emod_zz_vp, Eshell_xz_vp, Emod_xz_vp, Eshell_zx_vp, Emod_zx_vp] = deal(cell(1,length(indices)));
[lambda_shell_vp, lambda_mod_vp] = deal(cell(1,length(indices)));
Emod_xx_vp_mean.kmod = 0; Emod_xx_vp_mean.spec = 0;
Emod_zz_vp_mean.kmod = 0; Emod_zz_vp_mean.spec = 0;
Emod_xz_vp_mean.kmod = 0; Emod_xz_vp_mean.spec = 0;
Emod_zx_vp_mean.kmod = 0; Emod_zx_vp_mean.spec = 0;
c = 1;
for i = indices
    % Make symmetric along y
    Mvp2 = Mvp;
    Mvp2.Cxx = cat(3, flip(squeeze(Mvp.Cxx(:,:,:,i)), 3), squeeze(Mvp.Cxx(:,:,:,i)));
    Mvp2.Czz = cat(3, flip(squeeze(Mvp.Czz(:,:,:,i)), 3), squeeze(Mvp.Czz(:,:,:,i)));
    Mvp2.Cxz = cat(3, flip(squeeze(Mvp.Cxz(:,:,:,i)), 3), squeeze(Mvp.Cxz(:,:,:,i)));
    Mvp2.Czx = cat(3, flip(squeeze(Mvp.Czx(:,:,:,i)), 3), squeeze(Mvp.Czx(:,:,:,i)));
    Mvp2.dy = cat(2, -flip(Mvp.dy, 2), Mvp.dy);

    % Energy spectrum : xx
    [kx_vp, kz_vp, ky_vp, E_xx_vp] = compute3DFFT(Mvp2.dx, Mvp2.dz, Mvp2.dy, Mvp2.Cxx);
    [Eshell_xx_vp{c}, Emod_xx_vp{c}] = spectral_shells(kx_vp, kz_vp, ky_vp, E_xx_vp, N, kmax_vp, type);
    Emod_xx_vp_mean.kmod = ((c-1)*Emod_xx_vp_mean.kmod + Emod_xx_vp{c}.kmod)/c;
    Emod_xx_vp_mean.spec = ((c-1)*Emod_xx_vp_mean.spec + Emod_xx_vp{c}.spec)/c;

    % Energy spectrum : zz
    [~, ~, ~, E_zz_vp] = compute3DFFT(Mvp2.dx, Mvp2.dz, Mvp2.dy, Mvp2.Czz);
    [Eshell_zz_vp{c}, Emod_zz_vp{c}] = spectral_shells(kx_vp, kz_vp, ky_vp, E_zz_vp, N, kmax_vp, type);
    Emod_zz_vp_mean.kmod = ((c-1)*Emod_zz_vp_mean.kmod + Emod_zz_vp{c}.kmod)/c;
    Emod_zz_vp_mean.spec = ((c-1)*Emod_zz_vp_mean.spec + Emod_zz_vp{c}.spec)/c;
    
    % Energy spectrum : xz
    [~, ~, ~, E_xz_vp] = compute3DFFT(Mvp2.dx, Mvp2.dz, Mvp2.dy, Mvp2.Cxz(:,:,:));
    [Eshell_xz_vp{c}, Emod_xz_vp{c}] = spectral_shells(kx_vp, kz_vp, ky_vp, E_xz_vp, N, kmax_vp, type);
    Emod_xz_vp_mean.kmod = ((c-1)*Emod_xz_vp_mean.kmod + Emod_xz_vp{c}.kmod)/c;
    Emod_xz_vp_mean.spec = ((c-1)*Emod_xz_vp_mean.spec + Emod_xz_vp{c}.spec)/c;
    
    % Energy spectrum : zx
    [~, ~, ~, E_zx_vp] = compute3DFFT(Mvp2.dx, Mvp2.dz, Mvp2.dy, Mvp2.Czx(:,:,:));
    [Eshell_zx_vp{c}, Emod_zx_vp{c}] = spectral_shells(kx_vp, kz_vp, ky_vp, E_zx_vp, N, kmax_vp, type);
    Emod_zx_vp_mean.kmod = ((c-1)*Emod_zx_vp_mean.kmod + Emod_zx_vp{c}.kmod)/c;
    Emod_zx_vp_mean.spec = ((c-1)*Emod_zx_vp_mean.spec + Emod_zx_vp{c}.spec)/c;

    % Make xz-zx symmetric
    % E_xz_vp = 0.5*(E_xz_vp + conj(E_zx_vp));
    % E_zx_vp = conj(E_xz_vp);

    % Evaluate remaining elements of Energy spectrum matrix
    S = size(E_xx_vp);
    E_xy_vp = zeros(S); E_yx_vp = zeros(S); E_yz_vp = zeros(S); E_zy_vp = zeros(S); E_yy_vp = zeros(S);
    lambda_vp = cell(S);

    for j = 1:S(1)
        for k = 1:S(2)
            for l = 1:S(3)
                E_xy_vp(j,k,l) = - (kx_vp(j)/ky_vp(l))*E_xx_vp(j,k,l) - (kz_vp(k)/ky_vp(l))*E_xz_vp(j,k,l) ;
                E_yx_vp(j,k,l) = conj(E_xy_vp(j,k,l));

                E_zy_vp(j,k,l) = - (kx_vp(j)/ky_vp(l))*E_zx_vp(j,k,l) - (kz_vp(k)/ky_vp(l))*E_zz_vp(j,k,l) ;
                E_yz_vp(j,k,l) = conj(E_zy_vp(j,k,l));

                E_yy_vp(j,k,l) = - (kx_vp(j)/ky_vp(l))*E_yx_vp(j,k,l) - (kz_vp(k)/ky_vp(l))*E_yz_vp(j,k,l) ;

                E_matrix = [E_xx_vp(j,k,l) E_xy_vp(j,k,l) E_xz_vp(j,k,l) ; E_yx_vp(j,k,l) E_yy_vp(j,k,l) E_yz_vp(j,k,l) ; E_zx_vp(j,k,l) E_zy_vp(j,k,l) E_zz_vp(j,k,l)];
                eigval = eig(E_matrix);
                lambda_vp{j,k,l} = eigval;
            end
        end
    end

    [lambda_shell_vp{c}, lambda_mod_vp{c}] = analyze_lambda_ratio(kx_vp, kz_vp, ky_vp, lambda_vp, N, kmax_vp, type);

    c = c + 1;
end



style = {'--rs','--bo','--gd'};
figure(); 
subplot(2,1,1); hold on;
for i = 1:length(indices)
    plot(Emod_xx_vp{i}.kmod, Emod_xx_vp{i}.spec, style{i});
end
plot(Emod_xx_vp_mean.kmod, Emod_xx_vp_mean.spec, '--k', LineWidth=2);
legend('E_{xx}')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$E(k)$', 'Interpreter', 'latex');
yscale log;
subplot(2,1,2); hold on;
for i = 1:length(indices)
    plot(Emod_zz_vp{i}.kmod, Emod_zz_vp{i}.spec, style{i});
end
plot(Emod_zz_vp_mean.kmod, Emod_zz_vp_mean.spec, '--k', LineWidth=2);
legend('E_{zz}')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$E(k)$', 'Interpreter', 'latex');
yscale log;



style = {'--rs','--bo','--gd'};
figure(); 
subplot(2,1,1); hold on;
for i = 1:length(indices)
    plot(lambda_mod_vp{i}.kmod, lambda_mod_vp{i}.lambdaratio, style{i});
end
legend('\langle \lambda_{3}/\lambda_{2} \rangle')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$\langle \lambda_{3}/\lambda_{2} \rangle$', 'Interpreter', 'latex');
% yscale log;
subplot(2,1,2); hold on;
for i = 1:length(indices)
    plot(lambda_mod_vp{i}.kmod, lambda_mod_vp{i}.lambdamin, style{i});
end
legend('\langle \lambda_{2} \rangle')
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$\langle \lambda_{2} \rangle$', 'Interpreter', 'latex');
% yscale log;








M = Mvp;


for i = 1:length(M)
    [kx_vp, ky_vp, kz_vp, E_xx] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cxx);
    [Eshell_xx, Emod_xx] = spectral_shells(kx_vp, ky_vp, kz_vp, E_xx, N, [], 'log');

    [kx_yy, ky_yy, kz_yy, E_yy_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cyy);
    [Eshell_yy, Emod_yy] = spectral_shells(kx_yy, ky_yy, kz_yy, E_yy_vp, N, [], 'log');

    [kx_xy, ky_xy, kz_xy, E_xy_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cxy);
    [Eshell_xy, Emod_xy] = spectral_shells(kx_xy, ky_xy, kz_xy, E_xy_vp, N, [], 'log');

    [kx_yx, ky_yx, kz_yx, E_yx_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cyx);
    [Eshell_yx, Emod_yx] = spectral_shells(kx_yx, ky_yx, kz_yx, E_yx_vp, N, [], 'log');
end


figure(); hold on; 
plot(Emod_xx.kmod, Emod_xx.spec, '--rs'); 
plot(Emod_yy.kmod, Emod_yy.spec, '--bo');
xlabel('$k \ (mm^{-1})$', 'Interpreter', 'latex'); ylabel('$E(k)$', 'Interpreter', 'latex');
legend({'$E_{xx}$', '$E_{yy}$'},'Interpreter','latex')
xscale log; yscale log;


%%%%%%%%%%%% Horizontal Plane %%%%%%%%%%%%
load('C3D_SP_2021_0928.mat');




% saved .mat x is experiment x  .... Fig 1(a) in draft
% saved .mat y is experiment z
% saved .mat z is experiment y

[X, Y] = meshgrid(M(1).m.x, M(1).m.y);
figure(); contour(X, Y, squeeze(M(1).m.cxx(:,:,125))', 25);   


l0x = 10; l0y = 10; l0z = 10;

N = 128;


% Vertical Plane
for i = 1:length(M)
    [kx_vp, ky_vp, kz_vp, E_xx] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cxx);
    [Eshell_xx, Emod_xx] = spectral_shells(kx_vp, ky_vp, kz_vp, E_xx, N, [], 'log');

    [~, ~, ~, E_yy_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cyy);
    [Eshell_yy, Emod_yy] = spectral_shells(kx_vp, ky_vp, kz_vp, E_yy_vp, N, [], 'log');

    [~, ~, ~, E_xy_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cxy);
    [Eshell_xy, Emod_xy] = spectral_shells(kx_vp, ky_vp, kz_vp, E_xy_vp, N, [], 'log');

    [~, ~, ~, E_yx_vp] = compute3DFFT(M(i).m.x, M(i).m.y, M(i).m.z, M(i).m.cyx);
    [Eshell_yx, Emod_yx] = spectral_shells(kx_vp, ky_vp, kz_vp, E_yx_vp, N, [], 'log');

    figure(); plot(Emod_xx.kmod, Emod_xx.spec, 'rs'); xscale log; yscale log;


    % Evaluate remaining elements of Energy spectrum matrix
    S = size(E_xx);

    E_xz = zeros(S); E_zx = zeros(S); E_yz_vp = zeros(S); E_zy_vp = zeros(S); E_zz = zeros(S);
    lambda_vp = cell(S);

    for j = 1:S(1)
        for k = 1:S(1)
            for l = 1:S(1)
                E_xz(j,k,l) = - (kx_vp(j)/kz_vp(l))*E_xx(j,k,l) - (ky_vp(k)/kz_vp(l))*E_xy_vp(j,k,l) ;
                E_zx(j,k,l) = conj(E_xz(j,k,l));

                E_yz_vp(j,k,l) = - (kx_vp(j)/kz_vp(l))*E_yx_vp(j,k,l) - (ky_vp(k)/kz_vp(l))*E_yy_vp(j,k,l) ;
                E_zy_vp(j,k,l) = conj(E_yz_vp(j,k,l));

                E_zz(j,k,l) = - (kx_vp(j)/kz_vp(l))*E_zx(j,k,l) - (ky_vp(k)/kz_vp(l))*E_zy_vp(j,k,l) ;

                E_matrix = [E_xx(j,k,l) E_xy_vp(j,k,l) E_xz(j,k,l) ; E_yx_vp(j,k,l) E_yy_vp(j,k,l) E_yz_vp(j,k,l) ; E_zx(j,k,l) E_zy_vp(j,k,l) E_zz(j,k,l)];
                eigval = eig(E_matrix);
                lambda_vp{j,k,l} = eigval;
            end
        end
    end

    [lambda_shell_vp, lambda_mod_vp] = analyze_lambda_ratio(kx_vp, ky_vp, kz_vp, lambda_vp, N, [], 'log');

end



% Horizontal Plane



0.5*(2*pi/mean(diff(M(1).m.x)))
max(kx_vp)

2*pi/(max(abs(M(1).m.x))-min(abs(M(1).m.x)))
min(abs(kx_vp))