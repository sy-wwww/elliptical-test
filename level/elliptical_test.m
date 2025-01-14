function [T,sigma2,p_value] = elliptical_test(Xobs)
    [n,p] = size(Xobs);
    n1 = fix(n/2);
    n2 = n - n1;

    Xobs1 = Xobs(1:n1,:);
    Xobs2 = Xobs((n1+1):n,:);
        
    % data1 -- kappa1
    sigHat = (1/n1)*(Xobs1')*(Xobs1);
    kappa_p = transpose(mean(Xobs1.^4))./(diag(sigHat).^2)./3+ 2/n1;
    T1 = real(mean(kappa_p,'all'));

    % data2 -- kappa2
    sigHat = (1/n2)*(Xobs2')*(Xobs2);
    omega_n = trace(sigHat)^2;
    gamma_n = trace(sigHat^2)- omega_n/n2;

    l2norm = vecnorm(Xobs2').^2;
    sample_vec = sum((l2norm - mean(l2norm,'all')).^2,'all')/(n2-1);

    T2 = real((sample_vec +omega_n)/(omega_n + 2*gamma_n)); % fourth moment estimation
    
    
    T = (T1 - T2) * sqrt(n1*p);
    % all data -- sd
    sigHat = (1/n)*(Xobs')*(Xobs);
    b1 = trace(sigHat);
    b2 = trace(sigHat^2);
    b3 = trace(sigHat^3);
    b4 = trace(sigHat^4);
    
    var2 = (16*b4+8* (b2 - b1^2/n)^2)/ (b1^2+2*(b2 - b1^2/n))^2  *p;
    
    l2norm = vecnorm(Xobs').^2;
    sample_vec = sum((l2norm - mean(l2norm,'all')).^2,'all')/(n-1);
    a1h = 1 + real((sample_vec - 2*(b2 - b1^2/n))/(b1^2 + 2*(b2 - b1^2/n)));
    a2h = mean(vecnorm(Xobs').^6)./(b1^3+6*(b2 - b1^2/n)*b1+ 8 *b3);
    a3h = mean(vecnorm(Xobs').^8)/(b1^4+12*(b2 - b1^2/n)*b1^2+32*b1*b3 + 12 *(b2 - b1^2/n)^2 + 48 *b4);
    gammahat = a1h - 2*a2h;
    betahat = 1 - a3h;
    coe = 1 - betahat + gammahat; 
    if coe > log(p) * p^(-3/4)
        coe  = log(p) * p^(-3/4);
    elseif coe < -log(p) * p^(-3/4)
        coe  = -log(p) * p^(-3/4);     
    end
    
    A = real((1./sqrt(diag(sigHat)))'.*sigHat.*(1./sqrt(diag(sigHat))));
    var1 = 8/3*mean(A.^4,'all') * p * (1-betahat) + 8 * coe *mean(A.^2,'all')*p ;
    sigma2 = var1 + var2;
    
    p_value = normcdf(-abs(T/sqrt(sigma2)))*2;
end