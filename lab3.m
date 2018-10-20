% Laboratorio 3
display("Laboratorio 3")
printf("\n")
display("NOTA: EN ESTE LABORATORIO, LAS ABREVIACIONES CORRESPONDEN A:")
display("T1: MATRIZ DE TOEPLITZ DEL VECTOR CON N = 5")
display("T2: MATRIZ DE TOEPLITZ DEL VECTOR CON N = 25")
display("T3: MATRIZ DE TOEPLITZ DEL VECTOR DADO PARTE B")
display("Ta: MATRIZ DE VAMDER DEL VECTOR DADO PARTE A")
display("Tb: MATRIZ DE VAMDER DEL VECTOR DADO PARTE B")

% Antonella Requena 15-11196

% Ejercicio 1
% Los programas indicados en el enunciado
% estan ya implementados, se pueden encontrar
% en el mismo repositorio de este archivo

% Ejercicio 2
printf("\n")
display("Ejercicio 2")
% Resolucion de un sistema de ecuaciones Tx=b con T matriz de
% toeplitz

% Para n = 5
n = 11;   % es hasta donde llega el elemento 1/i del vector
V1 = [];  % se declara arreglo vacio

for i=1:1:n     % se llena el vector con los elementos
  p = 1/i;
  V1 = [V1,p];
end

T1 = toeplitz(V1);  % matriz de toeplitz

m1 = length(T1);
x1 = ones(m1,1);

b1 = T1*x1;

% Para n = 25
n = 51;
V2 = [];

for i=1:1:n
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
display("Autovalores de matriz de Toeplitz generada para n=5 T1")

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
display("Autovalores de matriz de Toeplitz generada para n=25 T2")

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
display("Autovalores de matriz de Toeplitz generada para V = −3, −2, −1, 0.01, 1, 2, 3 T3")

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

% 1era matriz de vander A

Ta = vander([0.5, 0.6, 0.7, 0.8, 0.9]);

x = ones(5,1);

b_a = Ta*x;

% 2da matriz de vander


Tb = vander([0.5, 0.6, 7, 8, 91013]);

b_b = Tb*x;


% Parte A

% Autovalores de Vander 1
display("Autovalores de matriz de Vander 1 generada")

eig(Ta)

% Determinante de Vander 1
printf("\n")
display("Determinante")
det(Ta)

% Numero de condicion de Vander 1
printf("\n")
display("Numero de condicion")
cond(Ta,inf)

% Autovalores de Vander 2
display("Autovalores de matriz de Vander 2 generada")

eig(Tb)

% Determinante de Vander 2
printf("\n")
display("Determinante")
det(Tb)

% Numero de condicion de Vander 2
printf("\n")
display("Numero de condicion")
cond(Tb,inf)


% Parte B

% Resolucion de sistemas de ecuaciones para Vander 1

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(Ta,b_a);

printf("\n")
display("Solucion del sistema para Vander 1 con Gauss sin pivoteo")
xa = sust_atras(A,b)

printf("\n")
display("Error relativo para Vander 1 con Gauss sin pivoteo")

ea= x-xa;
erelativo = norm(ea,inf)/norm(x,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(Ta,b_a);

printf("\n")
display("Solucion del sistema para Vander 1 con Gauss con pivoteo")
xa = sust_atras(A,b)

printf("\n")
display("Error relativo para Vander 1 con Gauss con pivoteo")

ea= x-xa;
erelativo = norm(ea,inf)/norm(x,inf)


% LU
printf("\n")
display("Solucion del sistema para Vander 1 con LU")
xa = sistemaLU(Ta,b_a)


printf("\n")
display("Error relativo para Vander 1 con LU")

ea= x-xa;
erelativo = norm(ea,inf)/norm(x,inf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Resolucion de sistemas de ecuaciones para T2

% Eliminacion de Gauss sin pivoteo

[A,b] = gauss(Tb,b_b);

printf("\n")
display("Solucion del sistema para Vander 2 con Gauss sin pivoteo")
xb = sust_atras(A,b)

printf("\n")
display("Error relativo para Vander 2 con Gauss sin pivoteo")

ea= x-xb;
erelativo = norm(ea,inf)/norm(x,inf)



% Eliminacion de Gauss con pivoteo
[A,b] = gauss_pivote(Tb,b_b);

printf("\n")
display("Solucion del sistema para Vander 2 con Gauss con pivoteo")
xb = sust_atras(A,b)

printf("\n")
display("Error relativo para Vander 2 con Gauss con pivoteo")

ea= x-xb;
erelativo = norm(ea,inf)/norm(x,inf)



% LU
printf("\n")
display("Solucion del sistema para Vander 2 con LU")
xb = sistemaLU(Tb,b_b)


printf("\n")
display("Error relativo para Vander 2 con LU")

ea= x-xb;
erelativo = norm(ea,inf)/norm(x,inf)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ejercicio 4
printf("\n")
display("Ejercicio 4")

% Calculamos inversa de T1
printf("\n")
display("Matriz T1 por su inversa")

inv = inversa(T1);        % Calculo de inversa con algoritmo eficiente
identidad =T1*inv        % Multiplicamos la matriz por su inversa

ea = identidad - eye(m1,m1); 

printf("\n")
display("Error absoluto de T1")

norm(ea,inf)            % error abs

% Calculamos inversa de T2
printf("\n")
display("Matriz T2 por su inversa")


inv = inversa(T2);
identidad =T2*inv

ea = identidad - eye(m2,m2);

printf("\n")
display("Error absoluto de T2")

norm(ea,inf)

% Calculamos inversa de T3
printf("\n")
display("Matriz T3 por su inversa")

inv = inversa(T3);
identidad =T3*inv

ea = identidad - eye(m3,m3);

printf("\n")
display("Error absoluto de T3")

norm(ea,inf)

% Calculamos inversa de Ta
printf("\n")
display("Matriz Ta por su inversa")

inv = inversa(Ta);
identidad =Ta*inv

ea = identidad - eye(5,5);

printf("\n")
display("Error absoluto de Ta")

norm(ea,inf)

% Calculamos inversa de Tb
printf("\n")
display("Matriz T3 por su inversa")


inv = inversa(Tb);
identidad =Tb*inv

ea = identidad - eye(5,5);
printf("\n")
display("Error absoluto de Tb")
norm(ea,inf)