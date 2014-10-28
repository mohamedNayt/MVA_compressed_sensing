%% Image Approximation with Orthogonal Bases %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
getd = @(p)path(p,path);
getd('../toolbox_signal/');
getd('../toolbox_general/');
%% Load image
n = 512;
f = rescale(load_image('lena',n));
imageplot(f);
%% Best M-terms Non-linear Approximation
%%
fF = fft2(f)/n;
clf;
imageplot(log(1e-5+abs(fftshift(fF))));
sortedfF = sort(abs(fF(:)),'descend');
%% M = N/100
M = floor(n^2/100);
T = sortedfF(M+1);
c = fF .* (abs(fF)>T);
fM = real(ifft2(c)*n);
SNR1 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));
%% M = N/20
M = floor(n^2/20);
T = sortedfF(M+1);
c = fF .* (abs(fF)>T);
fM = real(ifft2(c)*n);
SNR2 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));

%% Compression error as a function of M
semilogy(sortedfF);
err_fft = (sum(sortedfF.^2) - cumsum(sortedfF.^2))/sum(sortedfF.^2);
clf;
semilogy(err_fft);

%% Wavelet Approximation
%%
Jmin = 1;
options.h = compute_wavelet_filter('Daubechies',10);
fW = perform_wavortho_transf(f,Jmin,+1, options);
clf;
plot_wavelet(fW,Jmin);
title('Wavelet coefficients');
sortedfW = sort(abs(fW(:)),'descend');
%% M = N/100
M = floor(n^2/100);
T = sortedfW(M+1);
c = fW .* (abs(fW)>T);
fM = perform_wavortho_transf(c,Jmin,-1, options);
SNRW1 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));
%% M = N/20
M = floor(n^2/20);
T = sortedfW(M+1);
c = fW .* (abs(fW)>T);
fM = perform_wavortho_transf(c,Jmin,-1, options);
SNRW2 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));
%% fM convergence in Fourier vs Wavelet transform
err_wt = (sum(sortedfW.^2) - cumsum(sortedfW.^2))/sum(sortedfW.^2);
%semilogy(1:length(err_wt),err_wt,1:length(err_fft),err_fft);
semilogy([err_wt,err_fft]);
legend('Wavelets', 'Fourier');

%% Cosine Approximation
%%
fC = dct2(f);
clf;
imageplot(log(1e-5+abs(fC)));
sortedfC = sort(abs(fC(:)), 'descend');
%% M = N/100
M = floor(n^2/100);
T = sortedfC(M+1);
c = fC .* (abs(fC)>T);
fM = idct2(c);
SNRC1 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));
%% M = N/20
M = floor(n^2/20);
T = sortedfC(M+1);
c = fC .* (abs(fC)>T);
fM = idct2(c);
SNRC1 = -20*log10(norm(fM-f)/norm(f));
imageplot(clamp(fM));
%% fM convergence in DCT vs Fourier
err_dct = (sum(sortedfC.^2) - cumsum(sortedfC.^2))/sum(sortedfC.^2);
semilogy([err_dct,err_fft]);
legend('DCT', 'Fourier');

%% Local cosine approximation
%%
%% local cosine transform
w = 16;
fL = zeros(n,n);
for i = 1:n/w
    for j = 1:n/w
        iSel = w*(i-1)+1:w*i;
        jSel = w*(j-1)+1:w*j;
        P = f(iSel,jSel);
        fL(iSel,jSel) = dct2(P);
    end
end
clf;
imageplot(min(abs(fL),.005*w*w));
%% local cosine inverse transform
fL_inv = zeros(n,n);
for i = 1:n/w
    for j = 1:n/w
        iSel = w*(i-1)+1:w*i;
        jSel = w*(j-1)+1:w*j;
        P = fL(iSel,jSel);
        fL_inv(iSel,jSel) = idct2(P);
    end
end
clf;
imageplot(fL_inv);
%% M = N/100
sortedfL = sort(abs(fL(:)), 'descend');
M = floor(n^2/100);
T = sortedfL(M+1);
c = fL .* (abs(fL)>T);
fM = zeros(n,n);
for i = 1:n/w
    for j = 1:n/w
        iSel = w*(i-1)+1:w*i;
        jSel = w*(j-1)+1:w*j;
        P = c(iSel,jSel);
        fM(iSel,jSel) = idct2(P);
    end
end
clf;
imageplot(fM);
SNRC1 = -20*log10(norm(fM-f)/norm(f));
%% M = N/20
M = floor(n^2/20);
T = sortedfL(M+1);
c = fL .* (abs(fL)>T);
fM = zeros(n,n);
for i = 1:n/w
    for j = 1:n/w
        iSel = w*(i-1)+1:w*i;
        jSel = w*(j-1)+1:w*j;
        P = c(iSel,jSel);
        fM(iSel,jSel) = idct2(P);
    end
end
clf;
imageplot(fM);
SNRC2 = -20*log10(norm(fM-f)/norm(f));
%% fM convergence using local DCT
err_ldct = (sum(sortedfL.^2) - cumsum(abs(sortedfL).^2))/sum(sortedfL.^2);
xlimit = 7000;
semilogy([err_ldct(1:xlimit),err_dct(1:xlimit),err_wt(1:xlimit),err_fft(1:xlimit)]);
legend('Local DCT', 'DCT','Wavelets','Fourier');

%% Comparison of Wavelet Approximations of Several Images
%%
n = 512;
fList(:,:,1) = rescale( load_image('regular3',n) );
fList(:,:,2) = rescale( load_image('phantom',n) );
fList(:,:,3) = rescale( load_image('lena',n) );
fList(:,:,4) = rescale( load_image('mandrill',n) );
clf;
for i=1:4
    imageplot(fList(:,:,i),'', 2,2,i);
end
%% compute wavelet transform for each image
Jmin = 1;
options.h = compute_wavelet_filter('Daubechies',10);
err_wList = zeros(n^2,4);
for i=1:4
    fW = perform_wavortho_transf(fList(:,:,i),Jmin,+1, options);
    sortedfW2 = sort(abs(fW(:).^2),'descend');
    err_wList(:,i) = log10((sum(sortedfW2) - cumsum(sortedfW2))/sum(sortedfW2));
end
%% plot
xlimit = 7000;
x = log10(1:xlimit);
plot(x,[err_wList(1:xlimit,1),err_wList(1:xlimit,2),err_wList(1:xlimit,3),err_wList(1:xlimit,4)]);
legend('regular','phantom','lena','mandril');
