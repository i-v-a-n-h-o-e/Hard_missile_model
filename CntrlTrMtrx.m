function [Psi,Phi] = CntrlTrMtrx(A,B,Ts,varargin)
minargs = 3;        % Минимальное количество агрументов функции
maxargs = 4;        % Максимальное количество агрументов функции
narginchk(minargs, maxargs);    % Проверка количества входных параметров

if nargin == 3
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);  % Точность по умолчанию
elseif nargin == 4
    options = varargin{1};                          % Пользовательская точность
end

if size(A,1) == size(A,2)
    n = size(A,1);                              % Параметр размерности системы (количество строк матрицы А)
    if size(A,1) ~= size(B,1)
        error('ErrorTests:convertTest', ... 
              ['Матрица состояния и матрица управления/возмущений должны быть согласованы:'...
              '\nколичество строк в них должно совпадать']);
    end
    r = size(B,2);                              % Количество столбцов матрицы B
else
    error('Матрица состояния системы должна быть квадратной');
end
if (size(A,3) == 1)&&(size(B,3) == 1)
    % Вариант когда матрицы А и B стационарны (A(t) = const, B(t) = const)
    [Phi,Psi] = c2d(A,B,Ts);
elseif (size(A,3) == size(B,3))||((size(A,3) > 1) && (size(B,3) == 1))
    Phi = zeros(size(A));
    Psi = zeros([size(B,1),size(B,2),size(A,3)]);
    Phi(:,:,1) = eye(n);
    if size(B,3) > 1
        % Вариант когда A и В нестационарны (A(t) = F(t), B(t) = C(t))
        parfor k = 1:(size(A,3)-1)
            F = A(:,:,k:k+1);                   % F - два соседних слоя матрицы А: A(k) и A(k+1)
            C = B(:,:,k:k+1);                   % C - два соседних слоя матрицы B: B(k) и B(k+1)
            [~,zeta] = SolveEq2(F,C,Ts,options);
            Phi(:,:,k+1) = (reshape(zeta(end,1:n*n),n,n)).';
            Psi(:,:,k+1) = reshape(zeta(end,(n*n+1:n*n+n*r)),n,r);
        end
    else
        % Вариант когда A нестационарна, а B стационарна (A(t) = F(t), B(t) = const)
        C = zeros([size(B),2]);
        C(:,:,1) = B;
        C(:,:,2) = B;
        parfor k = 1:(size(A,3)-1)
            F = A(:,:,k:k+1);                   % F - два соседних слоя матрицы А: A(k) и A(k+1)
            [~,zeta] = SolveEq2(F,C,Ts,options);
            Phi(:,:,k+1) = (reshape(zeta(end,1:n*n),n,n)).';
            Psi(:,:,k+1) = reshape(zeta(end,n*n+1:n*r),n,r);
        end
    end
else
    error('Количество слоев матриц состояния и управления/возмущений не совпадают');
end

%% Функция решения дифференциального уравнения deta/dtau = -TET'(tau)*eta и интегрирование вектора ksi
function [T,ZETA] = SolveEq2(F,C,Ts,options)
% zeta - универсальный вектор, zeta = [eta; ksi], где eta = элементы
% матрицы Ф(t,tau), а ksi - 
n = size(F,1);                                  % Параметр размерности матрицы F
r = size(C,2);                                  % Количество столбцов матрицы C
IC = [reshape(eye(n),[],1);zeros(n*r,1)];       % Вектор начальных условий
[T,ZETA] = ode45(@(t,zeta) Eq2(t,zeta,Ts,F,C), [Ts,0], IC, options); % Здесь интегрирование производится в обратном времени

%% Функция для формирования дифференциального уравнения dzeta/dt = [-TET'(tau)*eta; -kappa]
function dzeta = Eq2(t,zeta,Ts,F,C)
% Eq2(t,y) - функция для записи дифференциального уравнения:
% dФ(t,tau)/dtau = -Ф(t,tau)*F(tau) в векторном виде deta/dtau = -TET(tau)*eta
n = size(F,1);                                  % Параметр размерности матрицы F
A = (F(:,:,2)-F(:,:,1))*t/Ts + F(:,:,1);        % Линейная интерполяция матрицы А на интервале [Ts*k,Ts*(k+1)]
B = (C(:,:,2)-C(:,:,1))*t/Ts + C(:,:,1);        % Линейная интерполяция матрицы B на интервале [Ts*k,Ts*(k+1)]
eta = zeta(1:n*n);
Phi = (reshape(eta,n,n)).';                     % Восстановление текущей матрицы Ф(t,tau)
kappa = reshape(Phi*B,[],1);                    % kappa - элементы матрицы подынтегрального выражения: Ф(t,tau)*C(tau)
TET = zeros(size(A,1)*size(A));
for k = 1:n                                     % Цикл формирования матрицы TET(tau)
    m = (k-1)*size(A,1) + 1;
    TET(m:k*size(A,1),m:k*size(A,2)) = A;
end
% Уравнение dzeta/dtau = [-TET'(tau)*eta; -kappa]
dzeta = [-TET.'*eta; -kappa];