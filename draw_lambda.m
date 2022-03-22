function draw_lambda(lambda, n)
    lambda_ = fliplr(lambda')
    N = [1:n];    
    epsilon = 2.22*10^(-16);

    plot(N,lambda_,'c+','Color','k')
    xlabel('j')
    ylabel('lambda_j(X)')
    axis([-inf, inf, epsilon * lambda_(1),lambda_(1)])
end