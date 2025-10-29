function dsp(x,c,intp_index,immin,immax)

if nargin==1
    c=0;
    intp_index=0;    
    immin=min(x(:));
    immax=max(x(:));
elseif nargin==2
    intp_index=0;    
    immin=min(x(:));
    immax=max(x(:));
elseif nargin==3
    immin=min(x(:));
    immax=max(x(:));    
end
if immin==immax
    sprintf('No display for constant image (intensity = %f)',immin)
    return
end
x=squeeze(x);
s=size(x);
figure
if s(2)==1
    ss=sqrt(s(1));
    x=reshape(x,ss,ss);
    if intp_index
        x=interp2(x,intp_index);
        immin=min(x(:));
        immax=max(x(:));
        imagesc(x,[immin immax]),axis equal,axis tight;
    else
        imagesc(x,[immin immax]),axis equal,axis tight;
    end
    if c==0
        colormap(gray),colorbar;
    else
        clinicalcolor,colorbar;
    end
elseif length(s)==2
    if s(1)>10*s(2)
        d=ceil(sqrt(s(2)));
        ss=sqrt(s(1));
        com_im=zeros(ss*d,ss*d);
        for n=1:s(2)
            [m,k]=ind2sub([d d],n);
            com_im((m-1)*ss+1:m*ss,(k-1)*ss+1:k*ss)=reshape(x(:,n),ss,ss);
        end
        if intp_index
            com_im=interp2(com_im,intp_index);
            immin=min(com_im(:));
            immax=max(com_im(:));
            imagesc(com_im,[immin immax]),axis equal,axis tight;
        else
            imagesc(com_im,[immin immax]),axis equal,axis tight;
        end
        if c==0
            colormap(gray),colorbar;
        else
            clinicalcolor,colorbar;
        end
    else
        if intp_index
            x=interp2(x,intp_index);
            immin=min(x(:));
            immax=max(x(:));
            imagesc(x,[immin immax]),axis equal,axis tight;
        else
            imagesc(x,[immin immax]),axis equal,axis tight;
        end
        if c==0
            colormap(gray),colorbar;
        else
            clinicalcolor,colorbar;
        end
    end
elseif length(s)==3
    d=ceil(sqrt(s(3)));
    com_im=zeros(s(1)*d,s(2)*d);
    for n=1:s(3)
        [m,k]=ind2sub([d d],n);
        com_im((m-1)*s(1)+1:m*s(1),(k-1)*s(2)+1:k*s(2))=x(:,:,n);
    end
    if intp_index
        com_im=interp2(com_im,intp_index);
        immin=min(com_im(:));
        immax=max(com_im(:));
        imagesc(com_im,[immin immax]),axis equal,axis tight;
    else
        imagesc(com_im,[immin immax]),axis equal,axis tight;
    end    
    if c==0
        colormap(gray),colorbar;
    else
        clinicalcolor,colorbar;
    end
end