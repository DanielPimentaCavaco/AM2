% NODE45 Método de ODE45 para resulução numérica de EDO/PVI

% INPUT:
%   f - función de EDO y'=f(t,y)
%   [a,b] - intervalo de valores da variável independente t
%   n - núnmero de subintervalos o interaciones del método
%   y0 - aproximación inicial y(a)=y0

%OUTPUT:
%   t - vector de intervalo [a,b] discretizado 
%   y - vector de las soluciones aproximadas de PVI en cada uno de los t(i)

function [t,y] = NODE45(f,a,b,n,y0)
    h = (b-a)/n;                % Tamanho de cada subintervalo (passo)
    t = a:h:b;                  % x
    [t,y] = ode45(f, t, y0);    % ODE45
    y = y';                     
end