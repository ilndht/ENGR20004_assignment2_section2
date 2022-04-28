close all 
clear all 
clc

%data for the beam at initial position
Before = readmatrix('Before_2.txt');
%data for the beam at P1 
P1 = readmatrix('P1_2.txt');
%data for the beam at P2 
P2 = readmatrix('P2_2.txt');
%data for the beam at P1+ P2
Both = readmatrix('Both_1.txt');

v = 0:1:292;

%Interpolated Before
xBefore = Before(:,1);
yBefore = Before(:,2);
newBefore = interp1(xBefore,yBefore,v);

% Interpolated P1
xP1 = P1(:,1);
yP1 = P1(:, 2);
newP1 = interp1(xP1,yP1,v,'spline');

%Interpolated P2
xP2 = P2(:,1);
yP2 = P2(:,2);
newP2 = interp1(xP2,yP2,v);

%Interpolated both
xBoth = Both(:,1);
yBoth = Both(:,2);
newBoth = interp1(xBoth,yBoth,v);

figure (1)
% the processed P1 
DefP1 = newP1 - newBefore;
plot(v,DefP1,'Linewidth',5)
title('Deflection curve of P1','FontSize',30)
xlim([0 292])
xlabel('x(mm)','FontSize',20);
ylabel('v(mm)','FontSize',20);
hold on

b = 20e-3 ; 
h = 1e-3 ; 

l = 0.073; 
F = 6.4 ; 
F2 = 2.2;
E = 69e9;
I = b*h^3/12 ; 
x = linspace(0,4*l,400); 
%analytical P1 
y = deflect(x,l,E,I,F)*1000;
plot(x*1000,y,'Linewidth',5);

plot(2*l*1000,0 ,'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 8)
legend('deflection from picture','analytical deflection','support C')

figure(2) 
%analytical P2 
y2 = deflect2(x,l,E,I,F2)*1000;
hold on

ylim([-30 30])
% the processed P2
DefP2 = newP2 - newBefore;
plot(v,DefP2,'Linewidth',5)
title('Deflection curve of P2','FontSize',30)
xlim([0 292])
plot(x*1000,y2,'Linewidth',5)
xlabel('x(mm)','FontSize',20);
ylabel('v(mm)','FontSize',20);
hold on
plot(2*l*1000,0 ,'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8)
legend('deflection from picture','analytical deflection','support C')


ylim([-30 30])

figure(3)
hold on
% the processed P1+P2
DefBoth = newBoth - newBefore;
plot(v,DefBoth,'-o' ,'Linewidth',5,'MarkerIndices', [73 73*3], 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'MarkerSize', 4)
title('Deflection curve of P1+P2 without superposition','FontSize',30)
xlim([0 292]);

xlabel('x(mm)','FontSize',20)
ylabel('v(mm)','FontSize',20)
DefBoth(1,73)
ylim([-30 30])
DefBoth(1, 73*3)
%analytical P1+P2 
plot(x*1000, y+y2,'Linewidth',5);
plot(2*l*1000,0 ,'ro','Linewidth',5, 'MarkerFaceColor', 'r', 'MarkerSize', 4)
legend('deflection from picture','analytical deflection','support C')

ylim([-30 30])

figure(4) 
%superposition P1+P2 
SuperP = DefP1 + DefP2;
plot(v, SuperP,'Linewidth',5)
hold on
% the processed picture of P1+P2 
plot(v,DefBoth,'Linewidth',5)
% analytical superposition 
plot(x*1000, y+y2,'Linewidth',5);
hold off
title('Deflection curve of P1 + P2 with superposition','FontSize',30)

xlim([0 292])
xlabel('x(mm)','FontSize',20)
ylabel('v(mm)','FontSize',20)

hold on 
ylim([-30 30])
%support C
plot(2*l*1000,0 ,'ro','Linewidth',5, 'MarkerFaceColor', 'r', 'MarkerSize', 4)
legend({'Superposed', 'picture','Analytical superposition','support C'})


% function for analytical P2 
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

% function for analytical P1
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
