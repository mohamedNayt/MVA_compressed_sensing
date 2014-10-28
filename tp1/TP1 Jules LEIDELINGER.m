getd = @(p)path(p,path);
getd('toolbox_signal/');
getd('toolbox_general/');

n = 512;
f = rescale( load_image('lena', n) );

clf;
imageplot(f);

fF = fft2(f)/n;
clf;
imageplot(log(1e-5 + abs(fftshift(fF))));

T = 20;
c = fF .* (abs(fF)>T);

fM = real(ifft2(c)*n);

imageplot(clamp(fM));

%exo1
N = n^2
M = floor(N/100)
%We find T by dichotomy, but it was not necessary to proceed like that...
Tmin = 0
Tmax = max(max(abs(fF)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fF)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fF)>Ttest));
end
c = fF .* (abs(fF)>Ttest);
fM = real(ifft2(c)*n);
imageplot(clamp(fM));

M = floor(N/20) 
Tmin = 0
Tmax = max(max(abs(fF)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fF)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fF)>Ttest));
end
c = fF .* (abs(fF)>Ttest);
fM = real(ifft2(c)*n);
figure;
imageplot(clamp(fM));

%exo2
cks = sort(abs(fF(:)));
cks = cks(end:-1:1);
semilogy(cks)

%exo3
ck2s = sort(abs(fF(:)).^2);
ck2s = ck2s(end:-1:1);
f2 = sum(ck2s);
err_fft = f2 - cumsum(ck2s);
semilogy(err_fft)

%Wavelet approximation

Jmin = 1;
options.h = compute_wavelet_filter('Daubechies',10);
fW = perform_wavortho_transf(f,Jmin,+1, options);

clf;
plot_wavelet(fW,Jmin);
title('Wavelet coefficients');

%exo4
N = n^2
M = floor(N/100)
%We find T by dichotomy, but it was not necessary to proceed like that...
Tmin = 0
Tmax = max(max(abs(fW)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fW)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fW)>Ttest));
end
c =fW .* (abs(fW)>Ttest);
fM = perform_wavortho_transf(c,Jmin,-1, options);
imageplot(fM);

M = floor(N/20) 
Tmin = 0
Tmax = max(max(abs(fW)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fW)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fW)>Ttest));
end
c = fW .* (abs(fW)>Ttest);
fM = perform_wavortho_transf(c,Jmin,-1, options);
figure;
imageplot(fM);

%exo5
cksW = sort(abs(fW(:)));
cksW = cksW(end:-1:1);
semilogy([cksW,cks])
legend('cksW','cks')

ck2sW = sort(abs(fW(:)).^2);
ck2sW = ck2sW(end:-1:1);
f2 = sum(ck2sW);
err_wav = f2 - cumsum(ck2sW);
semilogy([err_wav,err_fft])
legend('wavelets','fourier transform')

%Cosine approximation

fC = dct2(f);

clf;
imageplot(log(1e-5+abs(fC)));

%exo6
N = n^2
M = floor(N/100)
%We find T by dichotomy, but it was not necessary to proceed like that...
Tmin = 0
Tmax = max(max(abs(fC)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fC)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fC)>Ttest));
end
c = fC .* (abs(fC)>Ttest);
fM = real(idct2(c));
imageplot(clamp(fM));

M = floor(N/20) 
Tmin = 0
Tmax = max(max(abs(fC)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fC)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fC)>Ttest));
end
c = fC .* (abs(fC)>Ttest);
fM = real(idct2(c));
figure;
imageplot(clamp(fM));

%exo7
cksC = sort(abs(fC(:)));
cksC = cksC(end:-1:1);
semilogy([cksC,cksW,cks])
legend('cksC','cksW','cks')

ck2sC = sort(abs(fC(:)).^2);
ck2sC = ck2sC(end:-1:1);
f2 = sum(ck2sC);
err_dct = f2 - cumsum(ck2sC);
semilogy([err_dct,err_wav,err_fft])
legend('cosine','wavelets','fourier transform')

%Local Cosine Approximation

w = 16;

fL = zeros(n,n);

i = 5;
j = 7;

seli = (i-1)*w+1:i*w;
selj = (j-1)*w+1:j*w;
P = f(seli,selj);

fL(seli,selj) = dct2(P);

clf;
imageplot(P,'Patch',1,2,1);
imageplot(dct2(P-mean(P(:))),'DCT',1,2,2);

%exo8
for i = 1:n/w
    seli = (i-1)*w+1:i*w;
    for j = 1:n/w
        selj = (j-1)*w+1:j*w;
        P = f(seli,selj);
        fL(seli,selj) = dct2(P);
    end
end
clf;
imageplot(min(abs(fL),.005*w*w));

%exo9
for i = 1:n/w
    seli = (i-1)*w+1:i*w;
    for j = 1:n/w
        selj = (j-1)*w+1:j*w;
        P = fL(seli,selj);
        f1(seli,selj) = real(idct2(P));
    end
end
clf;
imageplot(f1);
error = norm(f(:)-f1(:),2)/norm(f(:),2);

%exo10
N = n^2
M = floor(N/100)
%We find T by dichotomy, but it was not necessary to proceed like that...
Tmin = 0
Tmax = max(max(abs(fL)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fL)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fL)>Ttest));
end
c = fL .* (abs(fL)>Ttest);
for i = 1:n/w
    seli = (i-1)*w+1:i*w;
    for j = 1:n/w
        selj = (j-1)*w+1:j*w;
        P = c(seli,selj);
        fM(seli,selj) = real(idct2(P));
    end
end
imageplot(clamp(fM));

M = floor(N/20) 
Tmin = 0
Tmax = max(max(abs(fL)))
Ttest = (Tmax + Tmin)/2
Mtest = sum(sum(abs(fL)>Ttest))
while Mtest ~= M
    if Mtest < M
        Tmax = Ttest;
    else
        Tmin = Ttest;
    end
    Ttest = (Tmax + Tmin)/2;
    Mtest = sum(sum(abs(fL)>Ttest));
end
c = fL .* (abs(fL)>Ttest);
for i = 1:n/w
    seli = (i-1)*w+1:i*w;
    for j = 1:n/w
        selj = (j-1)*w+1:j*w;
        P = c(seli,selj);
        fM(seli,selj) = real(idct2(P));
    end
end
imageplot(clamp(fM));

%exo11
cksL = sort(abs(fL(:)));
cksL = cksL(end:-1:1);
semilogy([cksL,cksC,cksW,cks])
legend('cksL','cksC','cksW','cks')

ck2sL = sort(abs(fL(:)).^2);
ck2sL = ck2sL(end:-1:1);
f2 = sum(ck2sL);
err_ldct = f2 - cumsum(ck2sL);
semilogy([err_ldct,err_dct,err_wav,err_fft])
legend('local cosine','cosine','wavelets','fourier transform')

%Comparison of Wavelet Approximations of Several Images
n = 512;
fList(:,:,1) = rescale( load_image('regular3',n) );
fList(:,:,2) = rescale( load_image('phantom',n) );
fList(:,:,3) = rescale( load_image('lena',n) );
fList(:,:,4) = rescale( load_image('mandrill',n) );

clf;
for i=1:4
    imageplot(fList(:,:,i),'', 2,2,i);
end

%exo12
fourier = [];
wavelet = [];
cosine = [];
localCos = [];
for index = [1 2 3 4]
    f = fList(:,:,index);
    %Fourier
    fF = fft2(f)/n;
    cfF = sort(abs(fF(:)));
    cfF = cfF(end:-1:1);
    %Wavelet
    Jmin = 1;
    options.h = compute_wavelet_filter('Daubechies',10);
    fW = perform_wavortho_transf(f,Jmin,+1, options);
    cfW = sort(abs(fW(:)));
    cfW = cfW(end:-1:1);
    %Cosine
    fC = dct2(f);
    cfC = sort(abs(fC(:)));
    cfC = cfC(end:-1:1);
    %Local Cosine
    w = 16;
    for i = 1:n/w
    seli = (i-1)*w+1:i*w;
        for j = 1:n/w
            selj = (j-1)*w+1:j*w;
            P = f(seli,selj);
            fL(seli,selj) = dct2(P);
        end
    end
    cfL = sort(abs(fL(:)));
    cfL = cfL(end:-1:1);
    row_fourier = [];
    row_wavelet = [];
    row_cosine = [];
    row_localCos = [];
    for M = unique(floor(logspace(0,5)))
        %Fourier
        T = cfF(M);
        c = fF .* (abs(fF)>=T);
        fM = real(ifft2(c)*n);
        err =norm(f(:)-fM(:))/norm(f(:));
        row_fourier = [row_fourier err];
        %Wavelet
        T = cfW(M);
        c = fW .* (abs(fW)>=T);
        fM = perform_wavortho_transf(c,Jmin,-1, options);
        err =norm(f(:)-fM(:))/norm(f(:));
        row_wavelet = [row_wavelet err];
        %Cosine
        T = cfC(M);
        c = fC .* (abs(fC)>=T);
        fM = real(idct2(c));
        err =norm(f(:)-fM(:))/norm(f(:));
        row_cosine = [row_cosine err];
        %Local Cosine
        T = cfL(M);
        c = fL .* (abs(fL)>=T);
        for i = 1:n/w
            seli = (i-1)*w+1:i*w;
            for j = 1:n/w
                selj = (j-1)*w+1:j*w;
                P = c(seli,selj);
                fM(seli,selj) = real(idct2(P));
            end
        end
        err =norm(f(:)-fM(:))/norm(f(:));
        row_localCos = [row_localCos err];
    end
    fourier = [fourier ; row_fourier];
    wavelet = [wavelet ; row_wavelet];
    cosine = [cosine ; row_cosine];
    localCos = [localCos ; row_localCos];
end

semilogy(transpose(fourier));
legend('regular','lena','phantom','mandrille')
title('Fourier')
figure;
semilogy(transpose(wavelet));
legend('regular','lena','phantom','mandrille')
title('Wavelet')
figure;
semilogy(transpose(cosine));
legend('regular','lena','phantom','mandrille')
title('Cosine')
figure;
semilogy(transpose(localCos));
legend('regular','lena','phantom','mandrille')
title('Local Cosine')
