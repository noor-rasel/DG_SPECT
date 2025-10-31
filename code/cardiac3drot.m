function c_im=cardiac3drot(x,axis_name)
%align original volume to LV axes (short, vertical and horizontal axes)
%NCAT: M<=>Y; N<=>X; S<=>Z
%Mingwu Jin, MIRC IIT
%May 09, 2006

if nargin<2
    axis_name='short';
end
% load 5dNCAT_orgIm x_roi
% x=mean(x_roi(:,:,:,1,1:2),5);
% x=x(21:52,14:45,:);
[M N S]=size(x);

%(1) rotate 55 degrees along Z-axis
for s=1:S
    c_im(:,:,s)=imrotate(x(:,:,s),55,'bilinear','crop');
end
%dsp(c_im): horizontal long axis (flipud needed!)
%dsp(permute(temp1,[2 3 1])); short
%dsp(permute(temp1,[3 1 2])); vertical long
if axis_name=='hlong'
    for s=1:S
        c_im(:,:,s)=flipud(c_im(:,:,s));
    end
    return
end
%(2) rotate along X-axis
c_im=permute(c_im,[3 1 2]);
for n=1:N
    c_im(:,:,n)=imrotate(c_im(:,:,n),-10,'bilinear','crop');
end
if axis_name=='vlong'
    for n=1:N
        c_im(:,:,n)=flipud(c_im(:,:,n));
    end
    return
end
%dsp(c_im2): vertical long axis (flipud needed!)
%(3) rotate along Y-axis
c_im=permute(c_im,[3 1 2]);
for m=1:M
    c_im(:,:,m)=imrotate(c_im(:,:,m),100,'bilinear','crop');
end
%dsp(c_im3): short axis