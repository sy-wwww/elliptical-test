function [n,p,tau,distr,Sigma,lambda_pop] = setting(code)
% sample size
n = 400;

% p/n
c1 = 0.5;
c2 = 1;
c3 = 1.5;
c4 = 3;
c5 = 5;
c6 = 10;

% distribution of \xi

distr1=@(n,p) sqrt(chi2rnd(p,n,1));
tau1 = 2;

tau2 = 3;
distr2=@(n,p) sqrt(beta_prime(n,p,tau2));


distr3 = @(n,p) sqrt((p+4).* betarnd(p/2,2,n,1)); 
tau3 = 0;


tau4 = 5;
distr4=@(n,p) sqrt(gamrnd(p/tau4,tau4,n,1)); 

distr5 = @(n,p) gamrnd(p,1,n,1) .* (1/sqrt(p+1));
tau5 = 4;



rho = 0.1;
k1 = 1/4;

switch code
	case 1
        c = c1; p=c*n; tau=tau1; distr=distr1;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5);
        Sigma = U*lambda_pop*U';
    case 2
        c = c2; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
	case 3
        c = c3; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 4
        c = c1; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 5
        c = c2; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 6
        c = c3; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 7
        c = c1; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 8
        c = c2; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 9
        c = c3; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 10
        c = c1; p=c*n; tau=tau2; distr=distr2;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 11
        c = c2; p=c*n; tau=tau2; distr=distr2;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 12
        c = c3; p=c*n; tau=tau2; distr=distr2;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 13
        c = c1; p=c*n; tau=tau2; distr=distr2;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 14
        c = c2; p=c*n; tau=tau2; distr=distr2;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;  
    case 15
        c = c3; p=c*n; tau=tau2; distr=distr2;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 16
        c = c1; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';      
    case 17
        c = c2; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 18
        c = c3; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 19
        c = c1; p=c*n; tau=tau3; distr=distr3;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 20
        c = c2; p=c*n; tau=tau3; distr=distr3;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 21
        c = c3; p=c*n; tau=tau3; distr=distr3;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 22
        c = c1; p=c*n; tau=tau3; distr=distr3;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 23
        c = c2; p=c*n; tau=tau3; distr=distr3;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;  
    case 24
        c = c3; p=c*n; tau=tau3; distr=distr3;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 25 
        c = c1; p=c*n; tau=tau3; distr=distr3;  
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 26 
        c = c2; p=c*n; tau=tau3; distr=distr3;  
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 27 
        c = c3; p=c*n; tau=tau3; distr=distr3; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 28
        c = c1; p=c*n; tau=tau4; distr=distr4;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 29
        c = c2; p=c*n; tau=tau4; distr=distr4;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 30
        c = c3; p=c*n; tau=tau4; distr=distr4;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 31
        c = c1; p=c*n; tau=tau4; distr=distr4;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 32
        c = c2; p=c*n; tau=tau4; distr=distr4;  
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 33
        c = c3; p=c*n; tau=tau4; distr=distr4; 
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 34 
        c = c1; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';  
    case 35 
        c = c2; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 36 
        c = c3; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 37
        c = c1; p=c*n; tau=tau5; distr=distr5;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 38
        c = c2; p=c*n; tau=tau5; distr=distr5;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 39
        c = c3; p=c*n; tau=tau5; distr=distr5;      
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 40
        c = c1; p=c*n; tau=tau5; distr=distr5;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1; 
    case 41
        c = c2; p=c*n; tau=tau5; distr=distr5;  
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1; 
    case 42
        c = c3; p=c*n; tau=tau5; distr=distr5;
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 43 
        c = c1; p=c*n; tau=tau5; distr=distr5; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 44 
        c = c2; p=c*n; tau=tau5; distr=distr5; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 45 
        c = c3; p=c*n; tau=tau5; distr=distr5; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';  
    % p/n =3,5
    case 46
        c = c4; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 47
        c = c4; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 48
        c = c4; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 49
        c = c5; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 50
        c = c5; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 51
        c = c5; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 52
        c = c4; p=c*n; tau=tau2; distr=distr2;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 53
        c = c4; p=c*n; tau=tau2; distr=distr2;     
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 54
        c = c4; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 55
        c = c5; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 56
        c = c5; p=c*n; tau=tau2; distr=distr2;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 57
        c = c5; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 58
        c = c4; p=c*n; tau=tau3; distr=distr3;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 59
        c = c4; p=c*n; tau=tau3; distr=distr3;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 60
        c = c4; p=c*n; tau=tau3; distr=distr3;  
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 61
        c = c5; p=c*n; tau=tau3; distr=distr3; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 62
        c = c5; p=c*n; tau=tau3; distr=distr3;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 63
        c = c5; p=c*n; tau=tau3; distr=distr3;  
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';       
    case 64
        c = c4; p=c*n; tau=tau4; distr=distr4;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 65
        c = c4; p=c*n; tau=tau4; distr=distr4;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 66
        c = c4; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 67
        c = c5; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 68
        c = c5; p=c*n; tau=tau4; distr=distr4;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 69
        c = c5; p=c*n; tau=tau4; distr=distr4;
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 70
        c = c4; p=c*n; tau=tau5; distr=distr5;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 71
        c = c4; p=c*n; tau=tau5; distr=distr5;     
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 72
        c = c4; p=c*n; tau=tau5; distr=distr5; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 73
        c = c5; p=c*n; tau=tau5; distr=distr5;
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 74
        c = c5; p=c*n; tau=tau5; distr=distr5;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 75
        c = c5; p=c*n; tau=tau5; distr=distr5;
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    % p/n =10
    case 76
        c = c6; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 77
        c = c6; p=c*n; tau=tau1; distr=distr1;   
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 78
        c = c6; p=c*n; tau=tau1; distr=distr1; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 79
        c = c6; p=c*n; tau=tau2; distr=distr2;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 80
        c = c6; p=c*n; tau=tau2; distr=distr2;     
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 81
        c = c6; p=c*n; tau=tau2; distr=distr2; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 82
        c = c6; p=c*n; tau=tau3; distr=distr3;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 83
        c = c6; p=c*n; tau=tau3; distr=distr3;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 84
        c = c6; p=c*n; tau=tau3; distr=distr3;  
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U'; 
    case 85
        c = c6; p=c*n; tau=tau4; distr=distr4;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 86
        c = c6; p=c*n; tau=tau4; distr=distr4;      
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 87
        c = c6; p=c*n; tau=tau4; distr=distr4; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
    case 88
        c = c6; p=c*n; tau=tau5; distr=distr5;   
        [U,R]=qr(randn(p));
        lambda_pop = eye(p);
        lambda_pop(1:5,1:5) = 5*eye(5); 
        Sigma = U*lambda_pop*U';
    case 89
        c = c6; p=c*n; tau=tau5; distr=distr5;     
        Sigma1 = zeros(p,p);
        for i = 1:p
           for j = 1:p
               Sigma1(i,j) = rho^abs(i-j);
           end
        end
        lambda_pop = diag(real(eig( Sigma1 )));
        Sigma = Sigma1;
    case 90
        c = c6; p=c*n; tau=tau5; distr=distr5; 
        [U,R]=qr(randn(p));
        lambda_pop = diag([1:p].^(-k1));
        Sigma = U*lambda_pop*U';
   % identity matrix     
    case 91   
        c = c1; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 92
        c = c2; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 93
        c = c3; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 94
        c = c4; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 95
        c = c5; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 96
        c = c6; p=c*n; tau=tau1; distr=distr1;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;   
    case 97   
        c = c1; p=c*n; tau=tau2; distr=distr2;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 98
        c = c2; p=c*n; tau=tau2; distr=distr2;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 99
        c = c3; p=c*n; tau=tau2; distr=distr2;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 100
        c = c4; p=c*n; tau=tau2; distr=distr2;      
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 101
        c = c5; p=c*n; tau=tau2; distr=distr2;     
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 102
        c = c6; p=c*n; tau=tau2; distr=distr2;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 103   
        c = c1; p=c*n; tau=tau3; distr=distr3;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 104
        c = c2; p=c*n; tau=tau3; distr=distr3;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 105
        c = c3; p=c*n; tau=tau3; distr=distr3;    
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 106
        c = c4; p=c*n; tau=tau3; distr=distr3;     
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 107
        c = c5; p=c*n; tau=tau3; distr=distr3;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 108
        c = c6; p=c*n; tau=tau3; distr=distr3;      
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 109   
        c = c1; p=c*n; tau=tau4; distr=distr4;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 110
        c = c2; p=c*n; tau=tau4; distr=distr4;    
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 111
        c = c3; p=c*n; tau=tau4; distr=distr4;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 112
        c = c4; p=c*n; tau=tau4; distr=distr4; 
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 113
        c = c5; p=c*n; tau=tau4; distr=distr4;     
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 114
        c = c6; p=c*n; tau=tau4; distr=distr4; 
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 115   
        c = c1; p=c*n; tau=tau5; distr=distr5;      
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 116
        c = c2; p=c*n; tau=tau5; distr=distr5;    
        lambda_pop = eye(p);
        Sigma = lambda_pop;
    case 117
        c = c3; p=c*n; tau=tau5; distr=distr5;    
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 118
        c = c4; p=c*n; tau=tau5; distr=distr5;   
        lambda_pop = eye(p);
        Sigma = lambda_pop; 
    case 119
        c = c5; p=c*n; tau=tau5; distr=distr5;   
        lambda_pop = eye(p);
        Sigma = lambda_pop;  
    case 120
        c = c6; p=c*n; tau=tau5; distr=distr5;     
        lambda_pop = eye(p);
        Sigma = lambda_pop;
end

end



