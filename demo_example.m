clear
close all
restoredefaultpath
addpath('./lib');

%% Demonstrate the code on a bi-Gaussian sphere distribution

a = 10;     % maximal radius
Nr = 25;    % bin number

% bin edges
r = linspace(0,a,Nr+1); r = r(:);       

% bin width
dr = mean(diff(r));

% bin centers
rc = r(1:end-1)/2 + r(2:end)/2;

% bi-Gaussian sphere radius sampling
rg1 = 6.5+1.2*randn(1e4,1);     % First Gaussian sampling
rg2 = 3.5+0.8*randn(2e4,1);     % Second Gaussian sampling
rg = [rg1;rg2];
rg(rg<0) = 0; rg(rg>a) = a;

% sphere radius histogram
Pv = histcounts(rg,r); Pv = Pv(:);
Pv = Pv/sum(Pv)/dr;

% cross-section radius histogram
Pa = Pv2Pa(Pv,rc,dr);

% Unfolding sphere radius histogram, gradient descent
xstarts = Pa;           % Initialize with cross-section radius histogram
lambda = 0;             % Regularization factor of L2 norm
weights = ones(Nr,1);   % Weights
% weights = 1./rc; weights = weights/sum(weights)*Nr;
options = optimoptions(@lsqnonlin,'Display','off','Algorithm','levenberg-marquardt');
Pv_fit = lsqnonlin(@(x)costfunction(Pa,x,rc,dr,lambda,weights),xstarts,[],[],options);

% Normalization
Pv_fit(Pv_fit<0) = 0;
Pv_fit = Pv_fit/sum(Pv_fit)/dr;

% Plot figure
figure; hold on;
hv = plot(rc,Pv,'-b','linewidth',1);            % Ground truth, sphere radius histogram
ha = plot(rc,Pa,'-r','linewidth',1);            % Input, cross-section radius histogram
hv_fit = plot(rc,Pv_fit,'--k','linewidth',1);   % Output, unfolded sphere radius histogram
legend([hv,ha,hv_fit],{'Ground truth, sphere radius histogram','Input, cross-section radius histogram','Output, unfolded sphere histogram'},'fontsize',20,'interpreter','latex');
box on; grid on;
set(gca,'fontsize',12);
xlabel('radius ($\mu$m)','interpreter','latex','fontsize',20);
ylabel('PDF ($\mu$m$^{-1}$)','interpreter','latex','fontsize',20);
xlim([0 a]); ylim([0 1/a*5]);






