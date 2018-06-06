clear all; clf;


%%Paramètres 
s = 20;%taille du cadre
h=0.2;

x0 = 0;
x1 = s;
y0 = 0;
y1 = s;
x = x0:h:x1;
y = y0:h:y1;
[X,Y] = meshgrid(x,y);
J1 = length(x);
J2 = length(y);
J=J1*J2;


%%Variables
D = 0.01;
Dc = 0.4;
alpha = 0.1;
lambda = 5;
c0 = 1;
cm = 40;
H = 10;
k = 1;
n0 = 1;


n = ones(J,1);
c = ones(J,1);

newn = ones(J,1);
newc = ones(J,1);






coinbasgauche = 1;
coinhautgauche = J1;
coinbasdroit = J1*(J2-1)+1;
coinhautdroit = J1*J2;
bordgauche = 2:J1-1;
borddroit = J1*(J2-2)+2 : J1*J2-1;
bordbas = J1+1:J1:J1*(J2-2)+1;
bordhaut = 2*J1-1 : J1 : J1*(J2-1)-1;


bord = [coinhautgauche, coinhautdroit, coinbasgauche, coinbasdroit, ...
     bordgauche, bordhaut, bordbas, borddroit];
interieur = setdiff(1:J, bord);


% interieur
L = sparse(interieur,interieur,-4,J,J); % matrice creuse, compacte en memoire
L = L + sparse(interieur,interieur+1,1,J,J);
L = L + sparse(interieur,interieur-1,1,J,J);
L = L + sparse(interieur,interieur+J2,1,J,J);
L = L + sparse(interieur,interieur-J2,1,J,J);




%%Conditions initiales
c(1:J)=1;
n(1:J)=0.9;


%%Conditions aux bords de Dirichlet
c(interieur)=1;
c(bord) = 1;
%c(borddroit) = 0.5;

n(bord)=1;
newn(bord)=1;
newc(bord)=1;


%Definition de la zone blessée


for i = 8/h:12/h
    n((8+i*s)/h+i:(12+i*s)/h+i)=0;
    c((8+i*s)/h+i:(12+i*s)/h+i)=0;
end

%%paramètres de simulation
tfinal = 10;
t0 = 0;
t = t0;
dt = 0.5*h^2/(4*max(D,Dc));
vect_t = t0:dt:tfinal;
vect_n = zeros(size(vect_t));
vect_c = zeros(size(vect_t));





beta = (c0^2+cm^2-2*h*c0*cm)/((c0-cm)^2);






%imagesc(N);

figure(1); clf;
%surf(X,Y,reshape(n,J,J),'EdgeColor','none');
surf(X,Y,reshape(c,J1,J2),'EdgeColor','none');
shading flat
view(3)
axis([0 20 0 20 0 6])
colorbar;
caxis([0,1]);
drawnow;
tk = 0;
pause

% BOUCLE PRINCIPALE
i=1;
while t < tfinal
    drawnow;
    
    
    %%%Activateur%%%
    
    %s = k*((2*cm*-H-beta)*c/(cm^2+c.^2))+beta;
    %f = (lambda*c0/n0)*((n0^2+alpha^2)/(n.^2+alpha^2))*n
    
    
    %%%Inhibiteur%%%
    
    %s = ((H-1).*c+H*c0)/(2*(H-1).*c+c0)*k
    %f = lambda*c0/n0.*n


   
     %
    
    newn =n + dt*  (D/(h^2)*L*n +  ((((H-1).*c+H*c0)./(2*(H-1).*c+c0)*k).*n).*(2-n/n0)      - k*n);
         
       
      
   
    newn(bord) = 1;
    
    
    
     %
    newc =c + dt* (D/(h^2)*L*c      + lambda*c0/n0.*n/3      - lambda * c);
      
       
    newc(bord)=1;
       
  
        
    n = newn;
    c = newc;
            
            
    %if tk > 1
        t
        %surf(X,Y,reshape(n,J1,J2),'EdgeColor','none');
        surf(X,Y,reshape(c,J1,J2),'EdgeColor','none');
        view(3)
        axis([0 20 0 20 0 10])
        colorbar;
        caxis([0,1]);
        drawnow;
        tk = 0;
    %end
    
    t = t + dt;
    %tk = tk + k;
    
    vect_n(i)=n(3*J1+10);
    vect_c(i)=c(3*J1+10);
    i=i+1;
    
    
    
    
    
    
end

plot(vect_t,vect_n,vect_t,vect_c);
legend('N','C');