%% **Image Deconvolution using Variational Method** %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
getd = @(p)path(p,path);
getd('../toolbox_signal/');
getd('../toolbox_general/');
%%
n = 256;
name = 'boat';
f0 = load_image(name);
f0 = rescale(crop(f0,n));
clf;
imageplot(f0);
%% Blur the image by convolution and noising
s = 3;
x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
h = exp( (-X.^2-Y.^2)/(2*s^2) );
h = h/sum(h(:));
hF = real(fft2(h));
clf;
imageplot(fftshift(h), 'Filter', 1,2,1);
imageplot(fftshift(hF), 'Fourier transform', 1,2,2);
%%
Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));
y0 = Phi(f0,h);
clf;
imageplot(f0, 'Image f0', 1,2,1);
imageplot(y0, 'Observation without noise', 1,2,2);
%%
sigma = .02;
y = y0 + randn(n)*sigma;
clf;
imageplot(y0, 'Observation without noise', 1,2,1);
imageplot(clamp(y), 'Observation with noise', 1,2,2);
%% Deconvulotion with L2 regularisation
%%
yF = fft2(y);
lambdas = 0:0.001:0.03;
SNR = zeros(length(lambdas),1);
for i=1:length(lambdas)
    lambda = lambdas(i);
    fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
    SNR(i) = snr(f0,fL2);
end
clf; plot(lambdas, SNR);
[m,ind] = max(SNR);
lambda = lambdas(ind);
%%
clf;
fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fL2), strcat(['L2 deconvolution, SNR=' num2str(snr(f0,fL2),3) 'dB']), 1,2,2);

%% Deconvolution by Sobolev Regularization.
%%
S = (X.^2 + Y.^2)*(2/n)^2;
lambdas = 0:0.001:0.25;
SNR = zeros(length(lambdas),1);
for i=1:length(lambdas)
    lambda = lambdas(i);
    fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
    SNR(i) = snr(f0,fSob);
end
clf; plot(lambdas, SNR);
[m,ind] = max(SNR);
lambda = lambdas(ind);
%%
fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fSob), strcat(['Sobolev deconvolution, SNR=' num2str(snr(f0,fSob),3) 'dB']), 1,2,2);

%% Deconvolution by Total Variation Regularization
%%
epsilon = 0.4*1e-2;
lambda = 0.06;
tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
niter = 600;
E = zeros(niter);

for i =1:niter
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );
    tv = sum(d(:));
    e = Phi(fTV,h)-y;
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
    E(i)  = 0.5*sum(e(:).^2) + lambda*tv;
end
clf;
plot(E)
imageplot(fTV)
%% Lambda selection
lambdas = 0.001:0.5*10^(-3):10^(-2);
SNR = zeros(length(lambdas));
for j=1:length(lambdas)
    lambda = lambdas(j);
    for i =1:niter
        Gr = grad(fTV);
        d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
        G = -div( Gr./repmat(d, [1 1 2])  );
        e = Phi(fTV,h)-y;
        fTV = fTV - tau*( Phi(e,h) + lambda*G);
    end
    SNR(j) = snr(f0,fTV);
end
