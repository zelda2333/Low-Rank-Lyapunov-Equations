% ADI solve Lyapunov
clc;
clear;

% N is the size of the matrix A and X
n = 4;
% Assign to the matrix A,B,X_0
h = 1 / (n+1);
epsilon = 0.1;
B = zeros([n,1]);
X_0 = zeros(n);
I = eye(n);
diag_0 = zeros([1,n]);
diag_0(1) = -1/h;
diag_1 = zeros([1,n-1]);

for i=2:n
    diag_0(i)= -2/h;
end
for i=1:n-1
    diag_1(i)= 1/h;
end

A_diag0 = diag(diag_0,0);
A_diag1 = diag(diag_1,1);
A_diag_1 = diag(diag_1,-1);

A = A_diag0 + A_diag1 + A_diag_1;
B(n) = 1/h;

% small->big
lambda = eig(-A); 
Im = sort(imag(lambda));
Re =sort(real(lambda));
% a = min(Re(lambda_i)); 
% b = max(Re(lambda_i));
% alpha = tan^-1 max|Im(lambda_i) / Re(lambda_i)|
a = Re(1);
b = Re(n);
% c=[1;2;3;4];%imag(lambda) 
Im_Re = sort(abs(imag(lambda) ./ real(lambda)));
alpha = tan(Im_Re(n))^(-1);
if Im==0
    dk = a / b;
else
    cos2_belta = 2 / (1 + 0.5 * (a / b + b / a));
    m = 2 * cos(alpha)^2 / cos2_belta -1;
    dk = (m + sqrt(m^2 - 1))^(-1);
end
k = sqrt(1-dk^2);
fK = @(x) sqrt(1 - k^2 * sin(x).^2).^(-1);
K = integral(fK, 0, 0.5*pi);
fv = @(x) sqrt(1 - dk^2 * sin(x).^2).^(-1);
v = integral(fv, 0, sin(sqrt(a / (b * dk)))^(-1));

J = ceil(K / (2 * v * pi) * log(4 / epsilon));
dn = zeros([1,J]);
p = zeros([1,J]);
P = zeros([J-1,n,n]);
z= zeros([n,1,J-1]);
% parameter P_J
for j = 1:J
    % dn(u,k) = sqrtm(1-k^2*sn(u,k)^2)
    dn(j) = sqrt(1 - k^2 * jacobiSN((2 * j - 1) * K / (2 * J), k)^2);
    p(j) = -sqrt(a * b / dk) * dn(j);    
end

for i = 1:J-1
   P(i,:,:) = (sqrt(-2 * p(i + 1)) / sqrt(-2 * p(i))) * (I - (p(i + 1) + p(i)) * (A + p(i + 1) * I)^(-1));
end

z(:,:,1) = sqrt(-2 * p(i)) * (A + p(1) * I)^(-1) * B;
Z_cfadi_j = z(:,:,1);
for j = 2:J
   z(:,:,j) = P(j-1) * z(j-1);
   if norm(z(:,:,j),2) > 10^(-8) || norm(z(:,:,j) / norm(z(:,:,j-1)),2) > 10^(-8) && j <= j
        Z_cfadi_ = cat(2, Z_cfadi_j, z(:,:,j));
        Z_cfadi_j = Z_cfadi_;
   end
end
Z_cfadi = Z_cfadi_
