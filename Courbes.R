rm(list =ls())

#Parametres :
a = 0.1
c0 = 10
n0 = 1
k= 10
l = 5
h = 2
cm = 20

B = (c0^2+cm^2-2*h*c0*cm)/((c0-cm)^2)

par(mfrow = c(2,1))

#dc/dt Activateur
curve(expr = l*c0*x/n0*((n0^2+a^2)/(x^2+a^2))-l*x, main = "c(n), k = , h = , alpha = , lambda = , cm = , n0 = , c0 = ", xlab = "n") 
#dn.dt Activateur
curve(expr = (x*k*(cm^2+c0^2)/(x^2+cm^2))*n0/(k*c0), main = "n(c), k = , h = , alpha = , lambda = , cm = , n0 = , c0 = ",xlab = "c",xlim = c(0,20)) 



#dc/dt inhibiteur
curve(expr = l*c0/n0*x-l*x, main = "c, k=10, lambda = 5, h=6, c0 = 1, n0 = 2",ylab = "f(n)-lambda*c", xlab = "n") 
#dn/dt inhibiteur
curve(expr = ((h-1)*x+h*c0)/(2*(h-1)*x+c0)*k, main = "n, k=10, lambda = 5, h=6, c0 = 1, n0 = 2",xlab = "c", ylab = "s(c)n(2-n/n0)", xlim = c(0,10)) 
