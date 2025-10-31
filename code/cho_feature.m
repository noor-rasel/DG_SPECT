function Imf=cho_feature(Fi,Im,loc)
%extract CHO features
%Feb. 9, 2007, Mingwu Jin

S=size(Fi);

Imf=zeros(S(3),1);
for i=1:S(3)
    iFi=ifft2(fftshift(Fi(:,:,i)));%,S(1),S(2));
    iFi=fftshift(iFi);
    iFi=wshift(2,iFi,loc);  % big ROI
    %iFi=abs(iFi);
    iFi=real(iFi);
    %stored ROI
    iFi=iFi(23:52,16:43);
    inter=iFi.*Im;
    Imf(i,1)=sum(inter(:));
end