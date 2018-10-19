% Laboratorio 3
% Antonella Requena 15-11196

% Ejercicio 1
% Los programas indicados en el enunciado
% estan ya implementados, se pueden encontrar
% en el mismo repositorio de este archivo

% Ejercicio 2
% Resolucion de un sistema de ecuaciones Tx=b con T matriz de
% toeplitz

% Para n = 5
n1 = 11;
V1 = [];

for i=1:1:n1
  p = 1/i;
  V1 = [V1,p];
end

T1 = toeplitz(V1);

m1 = length(T1);
x1 = ones(m1,1);

b1 = T1*x1;

% Para n = 25
n2 = 51;
V2 = [];

for i=1:1:n2
  p = 1/i;
  V2 = [V2,p];
end

T2 = toeplitz(V2);

m2 = length(T2);
x2 = ones(m2,1);

b2 = T2*x2;

% Para V = (−3, −2, −1, 0.01, 1, 2, 3)
V3 = [-3;-2;-1;0.01;1;2;3];

T3 = toeplitz(V3);

m3 = length(T3);
x3 = ones(m3,1);

b3 = T3*x3;

% Parte A

% Autovalores de T1
display("Autovalores de matriz de Toeplitz generada para n=5")

eig(T1)

% Determinante de T1
printf("\n")
display("Determinante")
det(T1)

% Numero de condicion de T1
printf("\n")
display("Numero de condicion")
cond(T1,inf)

% Autovalores de T2
display("Autovalores de matriz de Toeplitz generada para n=25")

eig(T2)

% Determinante de T2
printf("\n")
display("Determinante")
det(T2)

% Numero de condicion de T2
printf("\n")
display("Numero de condicion")
cond(T2,inf)

% Autovalores de T3
display("Autovalores de matriz de Toeplitz generada para V = −3, −2, −1, 0.01, 1, 2, 3")

eig(T3)

% Determinante de T3
printf("\n")
display("Determinante")
det(T3)

% Numero de condicion de T3
printf("\n")
display("Numero de condicion")
cond(T3,inf)

% Parte B

% Resolucion de sistemas de ecuaciones para T1

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(T1,b1);

printf("\n")
display("Solucion del sistema para T1 con Gauss sin pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T1 con Gauss sin pivoteo")

ea= x1-x;
erelativo = norm(ea,inf)/norm(x1,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(T1,b1);

printf("\n")
display("Solucion del sistema para T1 con Gauss con pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T1 con Gauss con pivoteo")

ea= x1-x;
erelativo = norm(ea,inf)/norm(x1,inf)



% LU
printf("\n")
display("Solucion del sistema para T1 con LU")
x = sistemaLU(T1,b1)


printf("\n")
display("Error relativo para T1 con LU")

ea= x1-x;
erelativo = norm(ea,inf)/norm(x1,inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolucion de sistemas de ecuaciones para T2

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(T2,b2);

printf("\n")
display("Solucion del sistema para T2 con Gauss sin pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T2 con Gauss sin pivoteo")

ea= x2-x;
erelativo = norm(ea,inf)/norm(x2,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(T2,b2);

printf("\n")
display("Solucion del sistema para T2 con Gauss con pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T2 con Gauss con pivoteo")

ea= x2-x;
erelativo = norm(ea,inf)/norm(x2,inf)



% LU
printf("\n")
display("Solucion del sistema para T2 con LU")
x = sistemaLU(T2,b2)


printf("\n")
display("Error relativo para T2 con LU")

ea= x2-x;
erelativo = norm(ea,inf)/norm(x2,inf)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolucion de sistemas de ecuaciones para T3

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(T3,b3);

printf("\n")
display("Solucion del sistema para T3 con Gauss sin pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T3 con Gauss sin pivoteo")

ea= x3-x;
erelativo = norm(ea,inf)/norm(x3,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(T3,b3);

printf("\n")
display("Solucion del sistema para T3 con Gauss con pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T3 con Gauss con pivoteo")

ea= x3-x;
erelativo = norm(ea,inf)/norm(x3,inf)



% LU

printf("\n")
display("Solucion del sistema para T3 con LU")
x = sistemaLU(T3,b3)


printf("\n")
display("Error relativo para T3 con LU")

ea= x3-x;
erelativo = norm(ea,inf)/norm(x3,inf)





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% Ejercicio 3
printf("\n")
display("Ejercicio 3")

% Resolucion de un sistema de ecuaciones Tx=b con T matriz de
% vander

% 1era matriz de vander

Ta = vander([0.5, 0.6, 0.7, 0.8, 0.9]);

xa = ones(5,1);

b1 = Ta*xa;

% 2da matriz de vander


Tb = vander([0.5, 0.6, 7, 8, 91013]);

xb = ones(5,1);

b2 = Tb*xb;


% Parte A

% Autovalores de T1
display("Autovalores de matriz de Vander 1 generada")

eig(Ta)

% Determinante de T1
printf("\n")
display("Determinante")
det(Ta)

% Numero de condicion de T1
printf("\n")
display("Numero de condicion")
cond(Ta,inf)

% Autovalores de T2
display("Autovalores de matriz de Vander 2 generada")

eig(Tb)

% Determinante de T2
printf("\n")
display("Determinante")
det(Tb)

% Numero de condicion de T2
printf("\n")
display("Numero de condicion")
cond(Tb,inf)


% Parte B

% Resolucion de sistemas de ecuaciones para T1

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(Ta,b1);

printf("\n")
display("Solucion del sistema para T1 con Gauss sin pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T1 con Gauss sin pivoteo")

ea= xa-x;
erelativo = norm(ea,inf)/norm(xa,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(Ta,b1);

printf("\n")
display("Solucion del sistema para T1 con Gauss con pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T1 con Gauss con pivoteo")

ea= xa-x;
erelativo = norm(ea,inf)/norm(xa,inf)



% LU
printf("\n")
display("Solucion del sistema para T1 con LU")
x = sistemaLU(Ta,b1)


printf("\n")
display("Error relativo para T1 con LU")

ea= xa-x;
erelativo = norm(ea,inf)/norm(xa,inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolucion de sistemas de ecuaciones para T2

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(Tb,b2);

printf("\n")
display("Solucion del sistema para T2 con Gauss sin pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T2 con Gauss sin pivoteo")

ea= xb-x;
erelativo = norm(ea,inf)/norm(xb,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(Tb,b2);

printf("\n")
display("Solucion del sistema para T2 con Gauss con pivoteo")
x = sust_atras(A,b)

printf("\n")
display("Error relativo para T2 con Gauss con pivoteo")

ea= xb-x;
erelativo = norm(ea,inf)/norm(xb,inf)



% LU
printf("\n")
display("Solucion del sistema para T2 con LU")
x = sistemaLU(Tb,b2)


printf("\n")
display("Error relativo para T2 con LU")

ea= xb-x;
erelativo = norm(ea,inf)/norm(xb,inf)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ejercicio 4
printf("\n")
display("Ejercicio 4")

% Calculamos inversa de T1
printf("\n")
display("Error abdoluto de T1")

inv = inversa_pro(T1);
identidad =T1*inv;

ea = identidad - eye(m1,m1);

norm(ea,inf)

% Calculamos inversa de T2
printf("\n")
display("Error abdoluto de T2")

inv = inversa_pro(T2);
identidad =T2*inv;

ea = identidad - eye(m2,m2);

norm(ea,inf)

% Calculamos inversa de T3
printf("\n")
display("Error abdoluto de T3")

inv = inversa_pro(T3);
identidad =T3*inv;

ea = identidad - eye(m3,m3);

norm(ea,inf)

% Calculamos inversa de Ta
printf("\n")
display("Error abdoluto de Ta")

inv = inversa_pro(Ta);
identidad =Ta*inv;

ea = identidad - eye(5,5);

norm(ea,inf)

% Calculamos inversa de Tb
printf("\n")
display("Error abdoluto de Tb")

inv = inversa_pro(Tb);
identidad =Tb*inv;

ea = identidad - eye(5,5);

norm(ea,inf)


















