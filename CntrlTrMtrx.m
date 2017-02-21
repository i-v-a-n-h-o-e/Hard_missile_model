function [Psi,Phi] = CntrlTrMtrx(A,B,Ts,varargin)
minargs = 3;        % ����������� ���������� ���������� �������
maxargs = 4;        % ������������ ���������� ���������� �������
narginchk(minargs, maxargs);    % �������� ���������� ������� ����������

if nargin == 3
    options = odeset('RelTol',1e-3,'AbsTol',1e-6);  % �������� �� ���������
elseif nargin == 4
    options = varargin{1};                          % ���������������� ��������
end

if size(A,1) == size(A,2)
    n = size(A,1);                              % �������� ����������� ������� (���������� ����� ������� �)
    if size(A,1) ~= size(B,1)
        error('ErrorTests:convertTest', ... 
              ['������� ��������� � ������� ����������/���������� ������ ���� �����������:'...
              '\n���������� ����� � ��� ������ ���������']);
    end
    r = size(B,2);                              % ���������� �������� ������� B
else
    error('������� ��������� ������� ������ ���� ����������');
end
if (size(A,3) == 1)&&(size(B,3) == 1)
    % ������� ����� ������� � � B ����������� (A(t) = const, B(t) = const)
    [Phi,Psi] = c2d(A,B,Ts);
elseif (size(A,3) == size(B,3))||((size(A,3) > 1) && (size(B,3) == 1))
    Phi = zeros(size(A));
    Psi = zeros([size(B,1),size(B,2),size(A,3)]);
    Phi(:,:,1) = eye(n);
    if size(B,3) > 1
        % ������� ����� A � � ������������� (A(t) = F(t), B(t) = C(t))
        parfor k = 1:(size(A,3)-1)
            F = A(:,:,k:k+1);                   % F - ��� �������� ���� ������� �: A(k) � A(k+1)
            C = B(:,:,k:k+1);                   % C - ��� �������� ���� ������� B: B(k) � B(k+1)
            [~,zeta] = SolveEq2(F,C,Ts,options);
            Phi(:,:,k+1) = (reshape(zeta(end,1:n*n),n,n)).';
            Psi(:,:,k+1) = reshape(zeta(end,(n*n+1:n*n+n*r)),n,r);
        end
    else
        % ������� ����� A �������������, � B ����������� (A(t) = F(t), B(t) = const)
        C = zeros([size(B),2]);
        C(:,:,1) = B;
        C(:,:,2) = B;
        parfor k = 1:(size(A,3)-1)
            F = A(:,:,k:k+1);                   % F - ��� �������� ���� ������� �: A(k) � A(k+1)
            [~,zeta] = SolveEq2(F,C,Ts,options);
            Phi(:,:,k+1) = (reshape(zeta(end,1:n*n),n,n)).';
            Psi(:,:,k+1) = reshape(zeta(end,n*n+1:n*r),n,r);
        end
    end
else
    error('���������� ����� ������ ��������� � ����������/���������� �� ���������');
end

%% ������� ������� ����������������� ��������� deta/dtau = -TET'(tau)*eta � �������������� ������� ksi
function [T,ZETA] = SolveEq2(F,C,Ts,options)
% zeta - ������������� ������, zeta = [eta; ksi], ��� eta = ��������
% ������� �(t,tau), � ksi - 
n = size(F,1);                                  % �������� ����������� ������� F
r = size(C,2);                                  % ���������� �������� ������� C
IC = [reshape(eye(n),[],1);zeros(n*r,1)];       % ������ ��������� �������
[T,ZETA] = ode45(@(t,zeta) Eq2(t,zeta,Ts,F,C), [Ts,0], IC, options); % ����� �������������� ������������ � �������� �������

%% ������� ��� ������������ ����������������� ��������� dzeta/dt = [-TET'(tau)*eta; -kappa]
function dzeta = Eq2(t,zeta,Ts,F,C)
% Eq2(t,y) - ������� ��� ������ ����������������� ���������:
% d�(t,tau)/dtau = -�(t,tau)*F(tau) � ��������� ���� deta/dtau = -TET(tau)*eta
n = size(F,1);                                  % �������� ����������� ������� F
A = (F(:,:,2)-F(:,:,1))*t/Ts + F(:,:,1);        % �������� ������������ ������� � �� ��������� [Ts*k,Ts*(k+1)]
B = (C(:,:,2)-C(:,:,1))*t/Ts + C(:,:,1);        % �������� ������������ ������� B �� ��������� [Ts*k,Ts*(k+1)]
eta = zeta(1:n*n);
Phi = (reshape(eta,n,n)).';                     % �������������� ������� ������� �(t,tau)
kappa = reshape(Phi*B,[],1);                    % kappa - �������� ������� ���������������� ���������: �(t,tau)*C(tau)
TET = zeros(size(A,1)*size(A));
for k = 1:n                                     % ���� ������������ ������� TET(tau)
    m = (k-1)*size(A,1) + 1;
    TET(m:k*size(A,1),m:k*size(A,2)) = A;
end
% ��������� dzeta/dtau = [-TET'(tau)*eta; -kappa]
dzeta = [-TET.'*eta; -kappa];