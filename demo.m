a = 3;
r = linspace(0,a,1e4);
rc = r(1:end-1)/2 + r(2:end)/2;
Pr = 1/a.*rc./sqrt(a^2-rc.^2);
dr = mean(diff(r));
sum(Pr.*dr)
figure; 

z = a*rand(1e6,1);
rs = sqrt(a^2-z.^2);
Nc = histcounts(rs,r);
Nc = Nc/sum(Nc)/dr;
hold on;
plot(rc,Nc,'-r');

plot(rc,Pr);
ylim([0 12])

%%
% close all

a = 10;
Nr = 25;
r = linspace(0,a,Nr+1); r = r(:);
dr = mean(diff(r));
rc = r(1:end-1)/2 + r(2:end)/2;
rg1 = 6.5+0.5*randn(1e4,1); 
rg2 = 2.5+0.8*randn(2e4,1);
rg = [rg1;rg2];
rg(rg<0) = 0; rg(rg>a) = a;
Nc = histcounts(rg,r); Nc = Nc(:);
Nc = Nc/sum(Nc)/dr;

tic;
Pa = Pv2Pa(Nc,rc,dr);
toc;

xstarts = Pa;
lambda = 0;
weights = ones(Nr,1);
% weights = 1./rc;
weights = weights/sum(weights)*Nr;
options = optimoptions(@lsqnonlin,'Display','off','Algorithm','levenberg-marquardt');%,'MaxIterations',2e3,'MaxFunctionEvaluations',2e3);
[xi,fxi] = lsqnonlin(@(x)costfunction(Pa,x,rc,dr,lambda,weights),xstarts,[],[],options);
xi(xi<0) = 0;
xi = xi/sum(xi)/dr;

figure; plot(rc,Nc,'-b'); hold on;
plot(rc,Pa,'-r');
plot(rc,xi,'--b');
Pa2 = Pv2Pa(xi,rc,dr);
plot(rc,Pa2,'--r');

%%
rg = sqrt(prod(cellsize,2));
a = 20;max(rg);
Nr = 20;
edges = linspace(0,a,Nr+1);
dr = mean(diff(edges));
rc = edges(1:end-1)/2 + edges(2:end)/2;
rc = rc(:);
Pa = histcounts(rg,edges);
Pa = Pa(:);
Pa = Pa/sum(Pa)/dr;



figure; 
bar(rc,Pa,0.3);
% plot(rc,Pa,'-r'); 

xstarts = Pa;
lambda = 0;
weights = ones(Nr,1);
weights = weights/sum(weights)*Nr;
options = optimoptions(@lsqnonlin,'Display','off','Algorithm','levenberg-marquardt');%,'MaxIterations',2e3,'MaxFunctionEvaluations',2e3);
[xi,fxi] = lsqnonlin(@(x)costfunction(Pa,x,rc,dr,lambda,weights),xstarts,[],[],options);
xi(xi<0) = 0;
xi = xi/sum(xi)/dr;
hold on; bar(rc+0.3,xi,0.3);
% hold on; plot(rc,xi,'-b');






