function empirical = elliptical_level(code,N)
    [n,p,tau,distr,Sigma,lambda_pop] = setting(code);
    pop_stat = zeros(N,3);
    for i=1:N
        xi = distr(n,p);
        X = randn(p,n);
        X = (transpose(xi)./vecnorm(X)).* X;
        Xobs = transpose(sqrtm(Sigma) * X);
        [T,sig,p_value] = elliptical_test(Xobs);
        pop_stat(i,:) = [T,sig,p_value];
    end
    empirical = mean(pop_stat(:,3)<0.05);
%     folder = fullfile( pwd);
%     save(strcat(folder,'/pop_stat_',int2str(code),'.txt'),'pop_stat','-ascii','-double');   
end