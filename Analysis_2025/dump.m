i = 1;
X = M(i).m.x; Y = M(i).m.y; Z = M(i).m.z; 

C = M(i).m.cxx;
Yplot = squeeze(C(:,61,125));
Yfit = log(Yplot);
figure(); plot(X, Yplot, '--rs'); hold on;
plot(X, exp(-1.0283)*exp(-X/3.68), '--k')
xlabel('$\Delta x$', 'Interpreter', 'latex');
ylabel('$\Gamma(\Delta x, \Delta y = 0, \Delta z = 0)$', 'Interpreter', 'latex');
yscale log

C = M(i).m.cyy;
Yplot = squeeze(C(56,:,125));
Yfit = log(Yplot);
figure(); plot(Y, Yplot, '--rs'); hold on;
plot(Y, exp(-0.5311)*exp(-Y/4.34), '--k')
xlabel('$\Delta y$', 'Interpreter', 'latex');
ylabel('$\Gamma(\Delta x = 0, \Delta y, \Delta z = 0)$', 'Interpreter', 'latex');
yscale log


