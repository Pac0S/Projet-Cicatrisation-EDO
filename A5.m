% Define and input parameters
% Note that $\lambda < 1/2$ to satisfy CFL condition.
% The code calculates an appropriate $dt$ given a lambda value.
D = 0.0005;
L = 1.;
a = 0.1; % alpha
c_m = 40;
h = 10;
lambda = 30;
D_C = 0.45;
b = (1+c_m*c_m-2*h*c_m)/((1-c_m)*(1-c_m)); % beta
N = 100;
tend = 30;
CFL = 0.2;
plotfreq = 1:1.5;
tplot = min(plotfreq,tend);
dx = L/(N+1);
dt = CFL*(dx*dx)/D_C;
b1 = D*dt/(dx*dx);
b2 = D_C*dt/(dx*dx);
x = dx*(0:1:N+1);
cold = zeros(N+2,1);
cnew = zeros(N+2,1);
nold = zeros(N+2,1);
nnew = zeros(N+2,1);
% Define initial conditions
for j=1:N+2
xx = x(j);
if xx < 1.0
nold(j) = 0;
cold(j) = 0;
else
nold(j) = 1.0;
cold(j) = 1.0;
end
end
% Plot horizontal line at .8 (wound declared healed at 80% of cell density).
nn = zeros(N+2,1);
for j=1:N+2
xx = x(j);
nn(j)=0.8;
end
plot(x,nn,'k:');
hold all;
% Define s(c) function.
s = @(x) (2*c_m*(h-b)*x)/(c_m*c_m+x*x)+b;
% Define f(n) function.
f = @(x) (x*(1+a*a))/(x*x+a*a);
% March forward in time using Forward Euler
t = dt;
tcount = 0;
while(t < tend)
% Set right-hand boundary condition (Dirichlet).
cnew(N+2) = 1.;
nnew(N+2) = 1.;
% March forward in time.
for j = N+1:-1:2
nnew(j) = nold(j) + b1*(nold(j-1)-2*nold(j)+nold(j+1))+dt*s(cold(j))*nold(j).*(2-nold(j))-nold(j)*dt;
cnew(j) = cold(j) + b2*(cold(j-1)-2*cold(j)+cold(j+1))+dt*lambda*(f(nnew(j))-cold(j));
end
% Set left-hand boundary condition (Neumann).
nnew(1) = nold(1) + b1*(-2*nold(1)+2*nold(2))+dt*s(cold(1))*nold(1)*(2-nold(1))-nold(1)*dt;
cnew(1) = cold(1) + b2*(-2*cold(1)+2*cold(2))+dt*lambda*(f(nold(1))-cold(1));
tcount = tcount+1;
% Plot solution curves.
if(t > tplot)
tplot = tplot+plotfreq;
plot(x,nnew,'g-')
plot(x,cnew,'b-')
end;
% Redefine for next iteration.
nold = nnew;
cold = cnew;
t = t + dt;
end;
% Plot all on same plot
axis([0 1 0 1.5]);
hold off;
name = strcat('Coupled Wound Healing System: CFL=',num2str(CFL),', dx=',num2str(dx),', dt=',num2str(dt));
xlabel('x')
ylabel('cell density n/ chemical concentration c')
title(name)