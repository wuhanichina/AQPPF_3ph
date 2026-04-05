function val = gk(mu, A, Sigma, k)
%GK   Compute the coefficient g^(k) for the quadratic form Q(X) = X'AX
%     where X ~ N(mu, Sigma).
%
% Usage:
%     val = gk(mu, A, Sigma, k)
%
% Inputs:
%     mu    : p-by-1 mean vector
%     A     : p-by-p matrix (usually symmetric)
%     Sigma : p-by-p covariance matrix
%     k     : non-negative integer (k = 0, 1, 2, ...)
%
% Output:
%     val   : The value of g^(k) = 2^k * k! * [ trace((A*Sigma)^(k+1)) 
%                                            + (k+1)* mu'*(A*Sigma)^k * A * mu ]
%
% References:
%   - Mathai, A. M., & Provost, S. B. (1992). Quadratic Forms in Random Variables. 
%   - Seber, G. A. F. (1977). Linear Regression Analysis.
%   - Kotz, S., Johnson, N. L., & Boyd, D. W. (1967). Ann. Math. Statist., 38, 832-848.

    % 简单的输入大小检查
    p = size(A, 1);
    assert(size(A,2) == p, 'Matrix A must be square.');
    assert(size(Sigma,1) == p && size(Sigma,2) == p, 'Sigma must be p-by-p.');
    assert(length(mu) == p, 'Length of mu must match dimensions of A and Sigma.');
    assert(k >= 0 && mod(k,1)==0, 'k must be a non-negative integer.');
    
    % 计算 (A*Sigma)^(k) 和 (A*Sigma)^(k+1)
    M = A*Sigma;
    Mk = M^k;          % (A*Sigma)^k
    Mk1 = M^(k+1);     % (A*Sigma)^(k+1)
    
    % 迹项 trace((A*Sigma)^(k+1))
    traceTerm = trace(Mk1);
    
    % 非中心项  (k+1)* mu'*(A*Sigma)^k * A * mu
    muTerm = (k+1)*(mu' * (Mk*A) * mu);
    
    % 带系数的拼合
    val = (2^k)*factorial(k) * (traceTerm + muTerm);
end