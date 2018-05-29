% Define and input parameters
% Note that $\lambda < 1/2$ to satisfy CFL condition.
% The code calculates an appropriate $dt$ given a lambda value.
D = 0.001;
a = 1.;
L = 1.;
N = 400;
tend = 400;
lambda = .25;
plotfreq = 1;
tplot = min(plotfreq,tend);
dx = L/(N+1);
dt = lambda*(dx*dx)/D;
b1 = D*dt/(dx*dx);
b2 = a*dt;
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
% March forward in time using Forward Euler
t = dt;
tcount = 0;
for t=1:12430
nnew(N+2) = 1;
for i=N+1:-1:2
nnew(i) = nold(i)+b1*(nold(i-1)-2*nold(i)+nold(i+1))+b2*nold(i).*(1-nold(i));
end
nnew(1) = nold(1)+b1*(-2*nold(1)+2*nold(2))+b2*nold(1).*(1-nold(1));
W(t)=find(nnew>.7999,1);
WW(t) = W(t)*(1/402);
tcount = tcount+1;
% Redefine for next iteration.
nold = nnew;
end;
% Plot time-wound radius graph
axis([0 1 0 1.2]);
WW(12431)=0;
WWW = [0:(1/12430):1];
plot(WWW, WW, 'k-')
name = strcat('Fishers Eq: D=',num2str(D),', lambda=',num2str(lambda),', dx=',num2str(dx),', dt=',num2str(dt))
xlabel('percentage of t')
ylabel('radius of wound')
title(name)