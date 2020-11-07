function mu_out = power_iters(A,sz)
%A = @(x)Aadj_3d(A3d(x));
%sz = [Ny, Nx, Nz];
bk = gpuArray(single(randn(sz)));
k = 0;
maxiter = 100;
mu = gpuArray(zeros(1));
resid = 1e10;
while resid > .0001 && k<=maxiter
    k = k+1;
    Abk =  A(bk);
    mu(k) = transpose(conj(bk(:)))*Abk(:)/(norm(bk(:))^2);
    if k >= 2
        resid = abs(mu(k) - mu(k-1));
    else
        resid = 1e10;
    end
    bknorm = norm(bk(:));
    bk = Abk/bknorm;
    fprintf('iter %i \t Eigenvalue %.6f \t residual %.6f \n',k,mu,resid)    
end
mu_out = mu(end);