getd = @(p)path(p,path); % scilab users must *not* execute this

getd('toolbox_signal/');
getd('toolbox_general/');

n = 256;
name = 'lena';
name = 'mri';
name = 'boat';
f0 = load_image(name);
f0 = rescale(crop(f0,n));

s = 3;

x = [0:n/2-1, -n/2:-1];
[Y,X] = meshgrid(x,x);
h = exp( (-X.^2-Y.^2)/(2*s^2) );
h = h/sum(h(:));

hF = real(fft2(h));

clf;
imageplot(fftshift(h), 'Filter', 1,2,1);
imageplot(fftshift(hF), 'Fourier transform', 1,2,2);

Phi = @(x,h)real(ifft2(fft2(x).*fft2(h)));

y0 = Phi(f0,h);

clf;
imageplot(f0, 'Image f0', 1,2,1);
imageplot(y0, 'Observation without noise', 1,2,2);

sigma = .02;

y = y0 + randn(n)*sigma;

clf;
imageplot(y0, 'Observation without noise', 1,2,1);
imageplot(clamp(y), 'Observation with noise', 1,2,2);

yF = fft2(y);

lambda = 0.02;

fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );

clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fL2), strcat(['L2 deconvolution, SNR=' num2str(snr(f0,fL2),3) 'dB']), 1,2,2);

%% exo1
errors = [];
lambdas =  0:0.0001:0.2;
for lambda = lambdas
    fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );
    errors = [errors snr(f0,fL2)];
end

plot(lambdas, errors);

[~, imax] = max(errors)
lambda = lambdas(imax)
fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );

clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fL2), strcat(['L2 deconvolution, SNR=' num2str(snr(f0,fL2),3) 'dB']), 1,2,2);

S = (X.^2 + Y.^2)*(2/n)^2;
lambda = 0.2;
fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
%%

%Deconvolution by Sobolev Regularization.

S = (X.^2 + Y.^2)*(2/n)^2;
lambda = 0.2;

fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );

clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fSob), strcat(['Sobolev deconvolution, SNR=' num2str(snr(f0,fSob),3) 'dB']), 1,2,2);

%% exo2
errors = [];
lambdas =  0:0.0001:0.2;
for lambda = lambdas
    fSob = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda*S) ) );
    errors = [errors snr(f0,fL2)];
end

plot(lambdas, errors);

[~, imax] = max(errors)
lambda = lambdas(imax)
fL2 = real( ifft2( yF .* hF ./ ( abs(hF).^2 + lambda) ) );

clf;
imageplot(y, strcat(['Observation, SNR=' num2str(snr(f0,y),3) 'dB']), 1,2,1);
imageplot(clamp(fSob), strcat(['Sobolev deconvolution, SNR=' num2str(snr(f0,fSob),3) 'dB']), 1,2,2);

%% Deconvolution
epsilon = 0.4*1e-2;
lambda = 0.06;
tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
niter = 600;

%% exo3

fTV = y;

tvs = [];

for iter = 1:niter
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );

    tvs = [tvs sum(d(:))];

    e = Phi(fTV,h)-y;
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
end

plot(tvs);

clf;
imageplot(clamp(fTV));


%% exo4

lambdas =  0:0.01:0.2;

errors = [];
for lambda = lambdas
    tau = 1.9 / ( 1 + lambda * 8 / epsilon);
    fTV = y;
    for iter = 1:niter
        Gr = grad(fTV);
        d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
        G = -div( Gr./repmat(d, [1 1 2])  );

        e = Phi(fTV,h)-y;
        fTV = fTV - tau*( Phi(e,h) + lambda*G);
    end
    errors = [errors snr(fTV,fL2)];
end

plot(errors);

[~, imax] = max(errors)
lambda = lambdas(imax)
tau = 1.9 / ( 1 + lambda * 8 / epsilon);
fTV = y;
for iter = 1:niter
    Gr = grad(fTV);
    d = sqrt( epsilon^2 + sum3(Gr.^2,3) );
    G = -div( Gr./repmat(d, [1 1 2])  );

    e = Phi(fTV,h)-y;
    fTV = fTV - tau*( Phi(e,h) + lambda*G);
end

clf;
imageplot(clamp(fSob), strcat(['Sobolev, SNR=' num2str(snr(f0,fSob),3) 'dB']), 1,2,1);
imageplot(clamp(fTV), strcat(['TV, SNR=' num2str(snr(f0,fTV),3) 'dB']), 1,2,2);
