% Some specific matrix as input data
clc;
clear;
% =========================================================================
% A CYCLIC LOW-RANK SMITH METHOD FOR LARGE SPARSE LYAPUNOV EQUATIONS
% Example 2.1 n = 400
% Example 6.1 n = 10,000
% =========================================================================
% N is the size of the matrix A and X
Input_1(400)
% Input_2(3600)
function Input_1(n)
    % Assign to the matrix A,B,ADI of epsilon,CFADI of tol
    epsilon = 2.22*10^(-16);
    tol = 10^(-8);
   
    h = 1 / (n + 1);
    B = zeros([n,1]);       
    diag_0 = zeros([1,n]);
    diag_0(1) = -1 / h;
    diag_1 = zeros([1, n-1]);

    for i=2:n
        diag_0(i)= -2 / h;
    end
    for i=1:n-1
        diag_1(i)= 1/h;
    end

    A_diag0 = diag(diag_0,0);
    A_diag1 = diag(diag_1,1);
    A_diag_1 = diag(diag_1,-1);

    A = A_diag0 + A_diag1 + A_diag_1;
    B(n) = 1 / h;
 
    X_ADI = ADI(A, B, n, epsilon);
    Z_cfadi = CFADI(A, B, n, epsilon, tol);
    X_CFADI = Z_cfadi * Z_cfadi'
    
end


% =========================================================================
% Example 6.2 n = 2n_0 = 3,000
% m<1
% Example 6.4 use reverse Cuthill–McKee algorithm(SYMRCM function),
% to reduce the bandwidth of matrix A, m=4
% =========================================================================

function Input_2(n0)
    % Assign to the matrix A,B,ADI of epsilon,CFADI of tol
    epsilon = 10^(-8);
    tol = 10^(-8);
   
    d = 1;
    k = 10;
    h = 1 / (n0 + 1);
    O = zeros(n0);
    I = eye(n0);
    B = zeros([n0,1]);       
    diag_0 = zeros([1,n0]);
    diag_0(1) = -k * h^(-2);
    diag_1 = zeros([1, n0-1]);

    for i=2:n0
        diag_0(i)= -2 * k * h^(-2);
    end
    for i=1:n0-1
        diag_1(i)= k * h^(-2);
    end

    A_diag0 = diag(diag_0,0);
    A_diag1 = diag(diag_1,1);
    A_diag_1 = diag(diag_1,-1);

    A_21 = A_diag0 + A_diag1 + A_diag_1;
    B(n0) = h^(-2);
    A = [O, I; A_21, -d * I];
 
%     % Example 6.4  
%     A = symrcm(A);
%     [i,j] = find(A);
%     % bandwidth of A
%     bwA = max(i-j) + 1;
    
    X_ADI = ADI(A, B, 2 * n0, epsilon)
    Z_cfadi = CFADI(A, B, 2 * n0, epsilon, tol);
    X_CFADI = Z_cfadi * Z_cfadi'
    
end




