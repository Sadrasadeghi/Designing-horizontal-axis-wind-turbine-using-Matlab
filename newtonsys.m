function [x] = newtonsys(F, JF, x0, tol, nmax)
%
% NEWTONSYS calcola una radice di un sistema non lineare
%
%       [x, R, niter] = newtonsys(F, JF, x0, tol, nmax)
%
% INPUT:
%   F      vettore che descrive le equazioni del sistema non lineare
%   JF     matrice contenente lo jacobiano di F
%   x0     vettore iniziale
%   nmax   numero di iterazioni massime
%
% OUTPUT:
%    ZERO     vettore radice di un sistema non lineare
%    RESIDUO  valore del residuo in ZERO
%    NITER    numero di iterazioni necessarie per calcolare ZERO
%    
% F e JF sono function sono definite tramite M-file, o come anonymous function


niter = 0;
err = tol + 1;
x = x0 ;
while (err >= tol && niter < nmax)
     Jx = JF(x); 
     R = F(x);  
     delta = - Jx\R ;
     x = x + delta ;
     err = norm ( delta );
     niter = niter + 1;
end
R = norm(F(x));
if ( niter == nmax && err > tol )
     fprintf([ 'Il metodo non converge nel massimo ' ,...
         ' numero di iterazioni . L ' ' ultima iterata \ n ' ,...
         ' calcolata ha residuo relativo pari a % e \ n '] , R );

end
return
