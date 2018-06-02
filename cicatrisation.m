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
J = length(x);
J2 = J*J;


n = ones(J2,1);
c = ones(J2,1);

newn = zeros(J2,1);
newc = zeros(J2,1);

%Definition de la zone blessée
for i = 8/h:12/h
    n((8+i*20)/h+i:(12+i*20)/h+i)=0;
    c((8+i*20)/h+i:(12+i*20)/h+i)=0;
end



%%paramètres de simulation
tfinal = 100000;
t0 = 0;
t = t0;
dt = 5;


%%Variables
D = 0.0005;
Dc = 0.45;
alpha = 0.1;
lambda = 30;
c0 = 1;
cm = 40;
H = 10;
k = 1;
n0 = 1;



beta = (c0^2+cm^2-2*h*c0*cm)/((c0-cm)^2);



% Condition periodiques
L = sparse(1:J2,1:J2,-2); % matrice creuse, compacte en memoire
coinhautgauche = 1;
coinbasgauche = J;
coinhautdroit = J*(J-1)+1;
coinbasdroit = J2;
bordgauche = 2:J-1;
bordhaut = J+1:J:J*(J-2)+1;
bordbas = 2*J:J:J*(J-1);
borddroit = J*(J-1)+2:J2-1;
bord = [coinhautgauche, coinhautdroit, coinbasgauche, coinbasdroit, ...
    bordgauche, bordhaut, bordbas, borddroit];
interieur = setdiff(1:J2, bord);

% interieur
L = L + sparse(interieur,interieur+1,1,J2,J2);
L = L + sparse(interieur,interieur-1,1,J2,J2);
L = L + sparse(interieur,interieur+J,1,J2,J2);
L = L + sparse(interieur,interieur-J,1,J2,J2);

% bords
L = L + sparse(bordhaut,bordhaut+1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J-1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J,1,J2,J2);
L = L + sparse(bordhaut,bordhaut-J,1,J2,J2);


L = L + sparse(bordgauche,bordgauche+1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche-1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J*(J-1),1,J2,J2);
% interieur
L = L + sparse(interieur,interieur+1,1,J2,J2);
L = L + sparse(interieur,interieur-1,1,J2,J2);
L = L + sparse(interieur,interieur+J,1,J2,J2);
L = L + sparse(interieur,interieur-J,1,J2,J2);

% bords
L = L + sparse(bordhaut,bordhaut+1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J-1,1,J2,J2);
L = L + sparse(bordhaut,bordhaut+J,1,J2,J2);
L = L + sparse(bordhaut,bordhaut-J,1,J2,J2);


L = L + sparse(bordgauche,bordgauche+1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche-1,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J,1,J2,J2);
L = L + sparse(bordgauche,bordgauche+J*(J-1),1,J2,J2);

L = L + sparse(bordbas,bordbas-(J-1),1,J2,J2);
L = L + sparse(bordbas,bordbas-1,1,J2,J2);
L = L + sparse(bordbas,bordbas+J,1,J2,J2);
L = L + sparse(bordbas,bordbas-J,1,J2,J2);

L = L + sparse(borddroit,borddroit+1,1,J2,J2);
L = L + sparse(borddroit,borddroit-1,1,J2,J2);
L = L + sparse(borddroit,borddroit-J*(J-1),1,J2,J2);
L = L + sparse(borddroit,borddroit-J,1,J2,J2);

% coins
L(coinhautgauche,coinhautgauche+1) = 1;
L(coinhautgauche,coinhautgauche+J-1) = 1;
L(coinhautgauche,coinhautgauche+J) = 1;
L(coinhautgauche,coinhautgauche+J*(J-1)) = 1;


L(coinbasgauche,coinbasgauche-(J-1)) = 1;
L(coinbasgauche,coinbasgauche-1) = 1;
L(coinbasgauche,coinbasgauche+J) = 1;
L(coinbasgauche,coinbasgauche+J*(J-1)) = 1;

L(coinhautdroit,coinhautdroit+1) = 1;
L(coinhautdroit,coinhautdroit+J-1) = 1;
L(coinhautdroit,coinhautdroit-J*(J-1)) = 1;
L(coinhautdroit,coinhautdroit-J) = 1;

L(coinbasdroit,coinbasdroit-(J-1)) = 1;
L(coinbasdroit,coinbasdroit-1) = 1;
L(coinbasdroit,coinbasdroit-J*(J-1)) = 1;
L(coinbasdroit,coinbasdroit-J) = 1;


%imagesc(N);

figure(1); clf;
%surf(X,Y,reshape(n,J,J),'EdgeColor','none');
surf(X,Y,reshape(c,J,J),'EdgeColor','none');
view(2)
drawnow;
tk = 0;
pause

% BOUCLE PRINCIPALE
while t < tfinal
    drawnow;
    
    %s = k*((2*cm*-H-beta)*c/(cm^2+c.^2))+beta;
    f = (lambda*c0/n0)*((n0^2+alpha^2)/(n.^2+alpha^2))*n


    newn(1:J,1) = 100;
    newn(J2-J:J2,1)=1;
   
    for i = 1:J-1
        newn(i*J+1,1)=1;
        newn(i*J,1)=1;
    end
    
    
    newc(1:J,1) = 100;
    newc(J2-J:J2,1)=1;
   
    for i = 1:J-1
        newc(i*J+1,1)=1;
        newc(i*J,1)=1;
    end
    
    
    newn = n + D*dt/(h^2)*L*n;% + ...
       %s*n - ...
       %k*n;
    
    newc = c + D*dt/(h^2)*L*c;% + ...
       %f*n - ...
       %lambda * c;
       
  
        
    n = newn;
    c = newc;
            
            
    %if tk > 1 
        t
        %surf(X,Y,reshape(n,J,J),'EdgeColor','none');
        surf(X,Y,reshape(c,J,J),'EdgeColor','none');
        view(2)
        drawnow;
        tk = 0;
    %end
    
    t = t + dt;
    %tk = tk + k;
    
    
end

