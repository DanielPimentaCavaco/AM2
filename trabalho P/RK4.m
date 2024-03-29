% NRK4 Método de Runge-Kutta de orden 2 para resolución numérica de EDO/PVI

% INPUT:
%   f - función de EDO y'=f(t,y)
%   [a,b] - intervalo de valores da variável independente t
%   n - núnmero de subintervalos o interaciones del método
%   y0 - aproximación inicial y(a)=y0

%OUTPUT:
%   t - vector de intervalo [a,b] discretizado 
%   y - vector de las soluciones aproximadas de PVI en cada uno de los t(i)

% 02/04/2024 - Alejandro Corral González .: a2023109587@isec.pt
% 02/04/2024 - Virginia Pérez Clemente .: a2023109270@isec.pt

function [t,y] = NRK4(f,a,b,n,y0)
h = (b-a)/n; % Tamaño de cada subintervalo (paso)

t = a:h:b;                                          % Asignación de memoria - vector de las abscisas
y = zeros(1, n+1);                                  % Asignación de memoria - vector de las ordenadas

y(1) = y0;                                          % El primer valor de y es siempre y0 (condición inicial del PVI)

for i=1:n                                           % El número de iteraciones será igual a n
    k1 = h*f(t(i), y(i));                           % Pendiente al inicio del intervalo
    k2 = h*f(t(i) + h/2, y(i) + 0.5*k1);            % Pendiente en el punto medio del intervalo
    k3 = h*f(t(i) + h/2, y(i) + 0.5*k2);            % Pendiente (nuevamente) en el punto medio del intervalo
    k4 = h*f(t(i+1), y(i) + k3);                    % Pendiente al final del intervalo
    
    y(i + 1) = y(i) + (k1 + 2*k2 + 2*k3 + k4)/6;    % Próximo valor aproximado de la solución del problema original
end
end