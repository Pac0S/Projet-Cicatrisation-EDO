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

N = ones(J,J);
C = ones(J,J);
N(18:32,22:28)=0;
C(18:32,22:28)=0;


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
tfinal = 200;
t0 = 0;
t = t0;
dt = 0.1;


%%Variables
D = 0.01;
Dc = 0.02;
lambda = 1;
cm = 2;
H = 10;
k = 0.1;
n0 = 1;
c0 = 0.5;
alpha = 0.1;

beta = (c0*c0+cm*cm-2*h*c0*cm)/((c0-cm)*(c0-cm));



% Condition periodiques
L = sparse(1:J2,1:J2,-4); % matrice creuse, compacte en memoire
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
surf(X,Y,reshape(n,J,J),'EdgeColor','none');
view(2)
drawnow;
tk = 0;
pause

% BOUCLE PRINCIPALE
while t < tfinal
    drawnow;
    %a =  k*((2*cm*(h-beta)*c(i)));
    %b = (cm*cm)+c(i)^2;
    
    %Ln = del2(n);
    %Lc = del2(c);
    %newn = n + D*Ln %+ k*((2*cm*(h-beta)*c)/(cm*cm+sum(c.*c))+beta) - k*n  ;
    %newc = c + Dc*Lc %+ lambda*c0*(n/n0)*((n0^2+alpha^2)/(sum(n.*n)+alpha^2)) - lambda * c;
    %newn = n + D*dt/(h^2)*L*n;
    
    newn = n + D*dt/(h^2)*L*n;% +k*c;%+ ...
       %k*((2*cm*(h-beta)*c)/(cm*cm)+c.^2) - ...
       %k*n;
    
    newc = c + Dc*dt/(h^2)*L*c;% + k*n;%+ ...
       %lambda*c0*n/n0*((n0^2+alpha^2)/(n.^2+alpha^2)) - ...
       %lambda * c;
        
    n = newn;
    c = newc;
            
            
    if tk > 1 
        t
        surf(X,Y,reshape(n,J,J),'EdgeColor','none');
        view(2)
        drawnow;
        tk = 0;
    end
    
    t = t + k;
    tk = tk + k;
    
    
end

