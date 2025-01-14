function rv = beta_prime(n,p,tau)
    xi_p = betarnd((1+p+tau)*p/tau, (1+p+2*tau)/tau,n,1);
    rv = xi_p./(1-xi_p);
end

