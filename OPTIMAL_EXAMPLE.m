% Some specific matrix as input data
clc;
clear;
% =========================================================================
% OPTIMAL ALTERNATING DIRECTION IMPLICIT PARAMETERS FOR NONSYMMETRIC SYSTEMS OF LINEAR EQUATIONS
% sigma <= 1, real; sigma >= 2,complex interval;1 < sigma < 2, Hermitian and skew-Hermitian parts
% sigma = 0, 1, 1.2, 1.4, 1.6, 1.8, 2, 3, 4, 5
% =========================================================================
% N is the size of the matrix A and X
Input_1(4,0.5)
function Input_1(n, sigma)
    % Assign to the matrix A,B,ADI of epsilon,CFADI of tol
    epsilon = 2.22*10^(-16);
    tol = 10^(-8);
   
    h = 1 / (n + 1);
    B = zeros([n,1]);       
    diag_0 = zeros([1,n]);
  
    diag_1 = zeros([1, n-1]);
    diag__1 = zeros([1, n-1]);

    for i=1:n-1
        diag_1(i)= -1 + sigma;
        diag__1(i)= -1 - sigma;
    end
    

    A_diag0 = diag(diag_0,0);
    A_diag1 = diag(diag_1,1);
    A_diag_1 = diag(diag__1,-1);

    A = 2 + A_diag0 + A_diag1 + A_diag_1;
    B(n) = 1 / h;
 
%     X_ADI = ADI(A, B, n, epsilon)
    Z_cfadi = CFADI(A, B, n, epsilon, tol);
    X_CFADI = Z_cfadi * Z_cfadi'
    
end