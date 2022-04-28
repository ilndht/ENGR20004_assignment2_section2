close all ; clear all ; clc
b = 20e-3 ; 
h = 1e-3 ; 

l = 0.073; 
F = 6.4 ; 
E = 69e9;
I = b*h^3/12 ; 
x = linspace(0,4*l,400); 

y = deflect(x,l,E,I,F);

plot(x,y)
hold on
plot(2*l,0 ,'bo')

function y = deflect(x,l,E,I,F)
y = zeros(size(x)) ; 
v1 =@(x) F/(E*I)*(-x.^3./12+x.*l^2./4); 
v2 = @(x) F/(2*E*I)*(x.^3./6 - l.*x.^2 + 1.5*l^2.*x-1/3*l^3);
v3 = @(x) -F*l^2.*x./(4*E*I) +F*l^3/(2*E*I) ; 
n = length(x) ; 
for i = 1: n 
    if x(i) <= l 
        y(i) = -v1(x(i)) ; 
    elseif x(i) <= 2*l && x(i)>= l 
        y(i) = -v2(x(i)) ; 
    else 
        y(i) = -v3(x(i)) ; 
    end 
end
end