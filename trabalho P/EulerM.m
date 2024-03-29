%NHEUN Método de Euler melhorado para resolução numérica de EDO/PVI

% Input:
%   f -> função do 2º menbro da ED
%   [a,b] -> extremos do intrevalo da variável independene t
%   n -> número de subintrevalos
%   y0 -> condição inicial (y=y0)

% Output:
%   y -> vector das soluções aproximantes
%   y(i+1) = y(i) + h*f(t(i),y(i)) , i=[0,n-1]

function [t,y] = EulerM(f,a,b,n,y0)
    h = (b-a)/n;
    t = a:h:b;          
    y = zeros(1,n+1);   
    y(1) = y0;
    
    for i=1:n                        
        k1 = f(t(i),y(i));           
        k2 = f(t(i+1), y(i) + k1*h); 
        k = 0.5*(k1+k2);             
        y(i+1)=y(i)+h*k;             
    end
end