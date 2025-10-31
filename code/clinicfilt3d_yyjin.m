% Erwan Gravier
% Research project
% 11/26/03
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implementation of the butterworth 3d
% filtering. Directly for use with the images
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Imfilt=clinicfilt3d(Image,order,cutoff)
%           Image: the 3D image to be filtered (64,64,64,k)
%                  size should be even
%           order: the order of the butterworth filter
%                  (should be integer eg 5)
%           cutoff: the cutoff frequency of the filter
%                   (should be less than 1, eg 0.2)
%           Imfilt: the 3D filtered image (64,64,64,k)
%                   # of counts kept
% Tested.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function    Imfilt=clinicfilt3d_yyjin(Image,order,cutoff)

% constants and variables
[ab1,ab2,ab3,ab4]=size(Image);
Tcounts=sum(Image(:));  % total number of counts
Imfreq=zeros(ab1,ab2,ab3);
Imfilt=zeros(ab1,ab2,ab3,ab4);

% creation of the filter
f=butter3d(ab1,cutoff,order);
ae=ab1/2;
f=circshift(f,[ae ae ae]);% Alignment on Frequency domain is necessary! -Mingwu
%fim=ifftn(f);  % put filter in image domain
%fim(fim<0)=0;  % remove  negative components (remove ringing effect) 
%f=ifftn(fim);  % put it back to fourrier domain, ready for use.

% filtering in the frequency domain
for i=1:ab4
    Imfreq=fftn(Image(:,:,:,i));
    Imfilt(:,:,:,i)=Imfreq.*f;
    Imfilt(:,:,:,i)=ifftn(Imfilt(:,:,:,i));
end

% postprocessing
Imfilt=real(Imfilt);  %take real part (avoid to get small imaginary parts...)

Imfilt(Imfilt<0)=0;     % remove negative values - YY

Imfilt=Imfilt/sum(Imfilt(:))*Tcounts;  % normalization, keep total number of counts




