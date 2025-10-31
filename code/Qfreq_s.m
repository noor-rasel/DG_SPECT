% Erwan Gravier
% Research project Jan 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function build 4 constant-Q
% frequency band.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Fi=Qfreq(S)  : 
%                S-> size of the frequency image of the filter (1 odd number)
%                Fi-> the 4 images of the filters in the frequency domain
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iFi=Qfreq_s(S,loc)

Fi=zeros(S,S,4);
alpha=2;
rhoc=S/(2*alpha^4);
mid=(S-1)/2+1;

for i=1:4
    rmin=floor(alpha^(i-1)*rhoc);
    rmax=floor(alpha^i*rhoc);
    for k=1:S
        for l=1:S
            r=sqrt((k-mid)^2+(l-mid)^2);
            if r>=rmin & r<rmax
                Fi(k,l,i)=1;
            end
        end
    end
end

for i=1:4
    iFi(:,:,i)=ifft2(fftshift(Fi(:,:,i)));%,S(1),S(2));
    iFi(:,:,i)=fftshift(iFi(:,:,i));
    iFi(:,:,i)=wshift(2,real(iFi(:,:,i)),loc);
end
%[-70 25]--> [199 104] = [129-loc(1) 129-loc(2)]