close all ; clear all ; clc
b = 20e-3 ; 
h = 1e-3 ; 

l = 0.073; 
F = 6.4 ; 
E = 69e9;
I = b*h^3/12 ; 
x = linspace(0,4*l,400); 

y2 = deflect2(x,l,E,I,F);

plot(x,y2)
hold on
plot(2*l,0 ,'bo')

function y = deflect2(x,l,E,I,F)
y = zeros(size(x)) ; 
v1 =@(x) F*x^3/(12*E*I)-F*l^2*x/(3*E*I); 
v2 = @(x) F/(E*I)*(-x^3/6+1.5*l*x^2-10/3*l^2*x+2*l^3);
v3 = @(x) 7/6*F*l^2*x/(E*I) -5/2*F*l^3/(E*I) ; 
n = length(x) ; 
for i = 1: n 
    if x(i) <= 2*l 
        y(i) = -v1(x(i)) ; 
    elseif x(i) <= 3*l && x(i)>= 2*l 
        y(i) = -v2(x(i)) ; 
    else 
        y(i) = -v3(x(i)) ; 
    end 
end


end