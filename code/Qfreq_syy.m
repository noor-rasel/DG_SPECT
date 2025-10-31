% Erwan Gravier
% Research project Jan 2003
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function build 4 constant-Q
% frequency band.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function iFi=Qfreq_syy(S,loc)


% Define band edges (normalized to Nyquist frequency)
bandedge = [2 4 8 16 32] / 70;
ch_num = length(bandedge) - 1;

Fi = zeros(S, S, ch_num);
mid = (S - 1)/2 + 1;

% Construct binary ring-shaped filters in frequency domain
for i = 1:ch_num
    rmin = floor(S * bandedge(i));
    rmax = floor(S * bandedge(i + 1));
    for k = 1:S
        for l = 1:S
            r = sqrt((k - mid)^2 + (l - mid)^2);
            if r >= rmin && r < rmax
                Fi(k, l, i) = 1;
            end
        end
    end
end


% Transform filters to spatial domain and shift to lesion location
for i = 1:ch_num
    iFi(:, :, i) = ifft2(fftshift(Fi(:, :, i)));
    iFi(:, :, i) = fftshift(iFi(:, :, i)); 
    iFi(:, :, i) = real(iFi(:, :, i));            
    iFi(:, :, i) = circshift(iFi(:, :, i), loc);  % Apply shift to loc
           
end  

%[-70 25]--> [199 104] = [129-loc(1) 129-loc(2)], 256*256 image coordinates