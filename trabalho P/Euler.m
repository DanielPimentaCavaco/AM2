function y = NEuler(f,a,b,n,y0)
% EULER Método de Euler para EDO/PVI
% 
% Y=Neuler(f,a,b,n,y0) Método para a resolição de um PVI
% y'= f(t,y) com t=[a,b] e y(a)=y0
% 
% Input:
%   f -> função do 2º menbro da ED
%   [a,b] -> extremos do intrevalo da variável independene t
%   n -> número de subintrevalos
%   y0 -> condição inicial (t=a logo y=y0)

% Output:
%   y -> vector das soluções aproximantes
%   y(i+1) = y(i) + h*f(t(i),y(i)) , i=[0,n-1]

h = (b-a)/n;
t=a:h:b;
y=zeros(1,n+1);
y(1) = y0;

    for i = 1:n
        y(i+1) = y(i)+h*f(t(i),y(i));
        t(i+1) = t(i)+h;
        end
end