function  [p, dp, chiqn, prob, corr]= fitpoli(nn,x,y,dy)
% funzione per il fit di un polinomio di grado n
% con la minimizzazione del chiquadro
% sintassi: [p, dp, chiqn, prob, corr]= fitpoli(n,x,y[,dy])
%    y = p(1)*x^(n) + p(2)*x^(n-1)+ ... p(n+1)
% I dati di ingresso sono i tre vettori x,y,dy e il grado del 
% polinomio
% Il dato dy puo' essere inserito sia come un vettore, 
% sia come un unico dato (in tal caso il vettore errori viene
% costruito con errori uguali per tutti i punti);
%
% Ritorna il vettore parametri, con il loro errore
% (corretto per il valore del chiquadro normalizzato),
% il chiquadro normalizzato, la probabilita' (valore alto ~ fit buono),
% e la matrice delle correlazioni
%
% ATTENZIONE: n=1 corrisponde alla retta, n=2 ad una parabola , ....
%             quindi la dimensione dei vettori P, DP, ..   e' N+1
%
% ADL 12.2005

n=nn+1;
if length(x)<n
    disp(' *** Numero di punti insufficiente! ***')
    disp(' Per determinare n coefficienti occorrono almeno n punti! ***')
    return
end
%----- controllo sulle dimensioni dei vettori------
if size(x) ~= size(y)
    disp(' *** I vettori x e y devono avere la stessa dimensione! ***')
    return
end
%-------------- controllo errori su y -------------
if nargin==3
    idy2=ones(size(x));
elseif size(dy)==1
    idy2=ones(size(x))/(dy.^2);
else
    idy2=1./(dy.^2);
end

%------- calcolo elementi per la matrice -----------
m=2*n-1;
sx=zeros(m,1);
for i=1:m
    sx(m-i+1) = sum(idy2 .* x.^(i-1));
end

%------- calcolo vettore termini noti -------------
sy=zeros(n,1);
for i=1:n
    sy(i)= sum(idy2 .* x.^(n-i).*y);
end

%------- calcolo elementi della matrice -----------
delta = zeros(n);
for i=1:n
    delta(:,i)=sx(i:i+n-1);
end

%------- calcolo inversa --------------------------
idelta = inv(delta);

%------- vettore dei parametri e errore -----------
p=idelta * sy;
dp= sqrt(diag(idelta));

%---- calcolo del chiquadro, della probabilita ----
res= y - polyval(p,x);
chiq=sum( idy2.*res.^2 );
nlib=length(x)-n;
chiqn=chiq/nlib;
prob=1-chi2cdf(chiq,nlib);
% 

%---- matrice correlazione ------------------------
corr = idelta;
for i=1:n
    for j=i:n
        corr(i,j)=corr(i,j)/(dp(i)*dp(j));
        corr(j,i)=corr(i,j);
    end
end

%---- correzione stima errori ---------------------
%dp=chiqn*dp;
return