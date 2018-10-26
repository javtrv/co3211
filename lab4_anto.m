% Laboratorio 4
% Antonella Requena 15-11196

% Ejercicio 1
display("Ejercicio 1")

A1 = [];
b1 = [];

for i=1:1:1000
    for j=1:1:1000
        if i==j
            A(i,j)=-5000;
        else
            A(i,j)= min(i,j);
        end
    end
end

b = [1:1:1000];

tic
A1,b1 = gauss(A,b);
x1 = sust_atras(A,b);
toc

tic
x2 = sistemaLUCrout(A,b);
toc

