clear all; clf;


%% Paramètres graphiques

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


%% Variables de simulation

D = 200;
Dc = 0.1;

H = 10;
k = 0.01;

cm = 40;
n0 = 1;
c0 = 1;

n = ones(J,1);
c = ones(J,1);

newn = ones(J,1);
newc = ones(J,1);

alpha = 0.1;
lambda = 10;

beta = (c0^2+cm^2-2*H*c0*cm)/((c0-cm)^2);




%% Laplacien discretisé

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



%{
%% Conditions aux bords de Dirichlet
c(bord) = 1;
n(bord)=1;
newn(bord)=1;
newc(bord)=1;
%}


%% Conditions initiales : zone blessée

% L = 8-2 = 6 
% l = 9-11 = 2

% dimensions de la blessures : 4*16 : (2:18, 8:12)

indices_blessure = ones(4/h,16/h).*(2/h*J1+8/h+1 : 2/h*J1 + 12/h)'; % Initialisation d'un vecteur aux dimentions de la blessures*
pas = 1:16/h; %Initialisation d'un vecteur de pas
indblessure = indices_blessure + J1* pas; %Indices des points de la blessure stockés dans une matrice

vectbles = indblessure(:); %Indices des points de la blessure stockés dans un vecteur


c(vectbles) = 0; %Zone blessée : concentrations nulles

%% Conditions initiales :Zone non blessée + condition de Dirichlet

nonbles = setdiff(1:J, vectbles); % Zone non blessée : concentrations à 1
c(nonbles) = 1;
n(nonbles)=1;


%{
%%  Bords de la blessure : conditions de Dirichlet

basgauche = 2/h*J1+8/h+1;
hautgauche = 2/h*J1 + 12/h;
basdroite = 18/h*J1+8/h+1;
hautdroite = 18/h*J1 + 12/h;


gauche =  basgauche+1:hautgauche-1;
droite =  basdroite+1:hautdroite-1;
bas = basgauche+J1 : J1 : basdroite - J1;
haut = hautgauche + J1 : J1 : hautdroite - J1;

matbordbles = [basgauche, gauche, hautgauche, bas, haut, basdroite, droite, hautdroite];
bordbles = matbordbles(:);

c(bordbles) = 1;
newc(bordbles) = 1;

%}

%% Paramètres de simulation
tfinal = 30;
t0 = 0;
t = t0;
dt = 0.5*h^2/(4*max(D,Dc));


%% Evolution de n et c dans la zone saine
vect_t = t0:dt:tfinal;
vect_n = zeros(size(vect_t));
vect_c = zeros(size(vect_t));



%% Affichage 

%imagesc(N);

figure(1); clf;
surf(X,Y,reshape(n,J1,J2),'EdgeColor','none');
surf(X,Y,reshape(c,J1,J2),'EdgeColor','none');
shading flat
view(2)
%view(3)
%axis([0 20 0 20 0 6])
colorbar;
caxis([0,1]);
drawnow;
tk = 0;
pause

%% BOUCLE PRINCIPALE
i=1;
while t < tfinal
    drawnow;
    
    
    %%%Activateur%%%
    
    %s = k*((2*cm*-H-beta)*c/(cm^2+c.^2))+beta;
    %f = (lambda*c0/n0)*((n0^2+alpha^2)/(n.^2+alpha^2))*n
    
    
    %%%Inhibiteur%%%
    
    %s = ((H-1).*c+H*c0)/(2*(H-1).*c+c0)*k
    %f = lambda*c0/n0.*n



    newn = n + dt*  (D/(h^2)*L*n +  ((((H-1).*c+H*c0)./(2*(H-1).*c+c0)*k).*n).*(2-n/n0)      - k*n);  
    newn(nonbles) = 1; %Dirichlet

    
    newc =c + dt* (D/(h^2)*L*c      + lambda*c0/n0.*n/3      - lambda * c);
    newc(nonbles) = 1; %Dirichlet
  
        
    n = newn;
    c = newc;
            
            
    %if tk > 1
        t
        surf(X,Y,reshape(n,J1,J2),'EdgeColor','none');
        surf(X,Y,reshape(c,J1,J2),'EdgeColor','none');
        view(2)
        %view(3)
        %axis([0 20 0 20 0 10])
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