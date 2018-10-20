% Laboratorio 3

% Antonella Requena

% Reporte de resultados

% Para cada caso los resultados fueron los siguientes: 

% Para la primera matriz de Toeplitz, a la cual llamamos T1 
% el error absoluto de la identidad que produce 
% T1 por su inversa es 2.6297e-09, siendo un numero
% bastante pequeño (casi 0) dado que la atriz es pequeña
% y sus elementos no tienen errores de representacion en el
% computador tan grandes (algunos son exactos)
% entonces al calcular la inversa, que se usa LU, el
% acarreo de errores es pequeño


% Para la segunda matriz de Toeplitz, a la cual llamamos T2,
% el error absoluto fue de 1.8008e+05, esto se debe a que,
% a diferencia de la matriz anterior, esta matriz es mucho 
% mas grande y sus elementos son números cada vez mas pequeños
% el calculo de la inversa acarrea más errores. Esto tambien
% puede reflejarse en que el numero de condición de la matriz es
% 6.1105e+18 lo que es un número grande. Es una matriz mal condicionada. 


% Para los casos de terca matriz de Toeplitz y primera matriz de 
% vander, los errores reportados fueron 
% Error absoluto de T3 : 4.6271e-12, error absoluto de Ta : 4.9880e-12
% ambos pequeños. Análisis análogo a la primera matriz.
% 

% Para la última matriz, la de vander B, el error reportado fue de 5.4913e+06
% esto puede explicarse dado a que la matriz de vander es una matriz mal condicionada
% su numero de condición es 1.8598e+21