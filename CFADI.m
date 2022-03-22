% CFADI solve Lyapunov

function Z_cfadi = CFADI(A, B, n, epsilon, tol)
    % small->big
    I = eye(n);
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
       if norm(z(:,:,j),2) > tol || norm(z(:,:,j) / norm(z(:,:,j-1)),2) > tol && j <= j
            Z_cfadi_ = cat(2, Z_cfadi_j, z(:,:,j));
            Z_cfadi_j = Z_cfadi_;
       end
    end
    Z_cfadi = Z_cfadi_;
end
