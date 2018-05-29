% Define and input parameters
% Note that $\lambda < 1/2$ to satisfy CFL condition.
% The code calculates an appropriate $dt$ given a lambda value.
D = 0.001;
s = 1.;
L = 1.;
N = 100;
tend = 30;
lambda = 0.3;
plotfreq = 1:3;
tplot = min(plotfreq,tend);
dx = L/(N+1);
dt = lambda*(dx*dx)/D;
b1 = D*dt/(dx*dx);
b2 = s*dt;
x = dx*(0:1:N+1);
nold = zeros(N+2,1);
nnew = zeros(N+2,1);
% Define initial conditions
for i=1:N+2
    xx = x(i);
    if xx < 1.0
nold(i) = 0.0;
    else
nold(i) = 1.0;
    end
end
% Plot initial conditions.
plot(x,nold,'r*')
hold all;
% Plot horizontal line at .8 (wound declared healed at 80% of cell density).
nn = zeros(N+2,1);
for i=1:N+2
    xx = x(i);
    nn(i)=0.8;
end
plot(x,nn,'k:');
hold all;
% March forward in time using Forward Euler
t = dt;
tcount = 0;
while(t < tend)
    nnew(N+2) = 1;
    for i=N+1:-1:2
        nnew(i) = nold(i)+b1*(nold(i-1)-2*nold(i)+nold(i+1))+b2*nold(i).*(1-nold(i));
    end
    nnew(1) = nold(1)+b1*(-2*nold(1)+2*nold(2))+b2*nold(1).*(1-nold(1));
    tcount = tcount+1;
% Plot solution curves at desired timesteps.
    if(t > tplot)
        tplot = tplot+plotfreq;
        plot(x,nnew,'b-')
    end;
% Redefine for next iteration.
    nold = nnew;
    t = t + dt;
end;
% Plot all on same plot
axis([0 1 0 1.2]);
hold off;
name = strcat('Fishers Eq: D=',num2str(D),', lambda=',num2str(lambda),', dx=',num2str(dx),', dt=',num2str(dt));
xlabel('x')
ylabel('cell density n')
title(name)