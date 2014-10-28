%% Inpainting using Sparse Regularization%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setup
getd = @(p)path(p,path);
getd('../toolbox_signal/');
getd('../toolbox_general/');

%% Sparse Regularization
%% Missing pixels and Inpainting
n = 128;
name = 'lena';
f0 = load_image(name);
f0 = rescale(crop(f0,n));
rho = .7;
Omega = ones(n,n);
sel = randperm(n^2);
Omega(sel(1:round(rho*n^2))) = 0;
Phi = @(f,Omega)f.*Omega;
clf;
imageplot(f0, 'Image f_0');
y = Phi(f0,Omega);
imageplot(clamp(y), 'Image y');
%% Soft Thresholding in a Basis
SoftThresh = @(x,T)x.*max( 0, 1-T./max(abs(x),1e-10) );
clf;
T = linspace(-1,1,1000);
plot( T, SoftThresh(T,.5) );
axis('equal');

Jmax = log2(n)-1;
Jmin = Jmax-3;

options.ti = 0; % use orthogonality.
Psi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
SoftThreshPsi = @(f,T)Psi(SoftThresh(PsiS(f),T));
clf;
imageplot( clamp(SoftThreshPsi(f0,.1)) );

%% Inpainting using Orthogonal Wavelet Sparsity
lambda = .03;
asVect = @(f)f(:);
E = @(f,lambda)0.5*norm(y - Omega.*f,'fro')^2 + lambda*sum(abs(asVect(PsiS(f))));

ProjC = @(f,Omega)(1-Omega).*f + y;
n_iter = 1000;
fSpars = y;
Es = zeros(n_iter);

for (i = 1:n_iter)
    Es(i) = E(fSpars,lambda);
    fSpars = ProjC(fSpars,Omega);
    fSpars = SoftThreshPsi( fSpars, lambda );
end
clf;
plot(Es);
imageplot(clamp(fSpars));

%% Inpainting using Translation Invariant Wavelet Sparsity
J = Jmax-Jmin+1;
u = [4^(-J) 4.^(-floor(J+2/3:-1/3:1)) ];
U = repmat( reshape(u,[1 1 length(u)]), [n n 1] );
lambda = .01;
options.ti = 1; % use translation invariance
Xi = @(a)perform_wavelet_transf(a, Jmin, -1,options);
PsiS = @(f)perform_wavelet_transf(f, Jmin, +1,options);
Psi = @(a)Xi(a./U);
E2 = @(a, lambda)0.5*norm(y-Phi(Psi(a), Omega), 'fro') + lambda*sum(abs(asVect(a)));
%init
n_iter2 = 1000;
Es2 = zeros(n_iter2);
a = U.*PsiS(fSpars);
%loop
for i = 1:n_iter
   Es2(i) = E2(a, lambda);
   tau = 1.9*min(u);
   fTI = Psi(a);
   a = a + tau*PsiS( Phi( y-Phi(fTI,Omega),Omega ) );
   a = SoftThresh( a, lambda*tau ); 
end
plot(Es2);
clf;
imageplot(clamp(Psi(a)));
