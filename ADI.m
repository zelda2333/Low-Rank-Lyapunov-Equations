% ADI solve Lyapunov

function X_ADI = ADI(A, B, n, epsilon)
    X_0 = zeros(n); 
    I = eye(n);
    % small->big
    lambda = eig(-A);
    % draw_lambda
%     draw_lambda(lambda, n)
    
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
    P = zeros([1,J]);

    % parameter P_J
    for j=1:J
        % dn(u,k) = sqrtm(1-k^2*sn(u,k)^2)
        dn(j) = sqrt(1 - k^2 * jacobiSN((2 * j - 1) * K / (2 * J), k)^2);
        P(j) = -sqrt(a * b / dk) * dn(j);
        X_j2 = -(A + P(j) * I)^(-1) * (B * B' + X_0 * (A' - P(j) * I));
        X_j = -(A + P(j) * I)^(-1) * (B * B' + X_j2' * (A' - P(j) * I));

        X_0 = X_j;
    end

    X_ADI = X_0;
end
