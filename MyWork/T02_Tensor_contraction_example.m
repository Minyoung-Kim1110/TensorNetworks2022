d_a = 101; % d_alpha
d_b = 102; % d_beta
d_c = 103; % d_gamma
d_d = 104; % d_delta
d_m = 105; % d_mu

A = rand(d_c,d_d);     % tensor A(gamma,delta)
B = rand(d_a,d_m,d_c); % tensor B(alpha,mu,gamma)
C = rand(d_b,d_m,d_d); % tensor C(beta,mu,delta)


tobj = tic2;

% A, C contract 
C1 = permute(C,[3,1,2]);
C1 = reshape(C1, [d_d, d_b * d_m]);
AC = A * C1 ; % tensor AC (gamma, beta * mu)
AC = reshape(AC, [d_c, d_b, d_m]);

% B contract with AC 
AC = permute(AC, [3,1,2]);
AC = reshape(AC, [d_m* d_c, d_b]);
B1 = reshape(B, [d_a, d_m * d_c]);
BAC = B1 * AC ; 

% compare real time and cpu time 
toc2(tobj,'-v');
% real time << cpu time due to multi threading 

% check
% whos BAC;



