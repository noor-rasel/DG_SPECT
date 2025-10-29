function cor = template_matching(template,target_image,mm,x,y,z,gg)

[x_template,y_template,z_template] = size(template);
[x_image,y_image,z_image] = size(target_image);

x_pad = (x_template-1)/2;
y_pad = (y_template-1)/2;
z_pad = (z_template-1)/2;

target_pad = padarray(target_image,[x_pad,y_pad,z_pad]);

if gg<4
    x_start = x_pad+x-4;
    y_start = y_pad+y-4;
    z_start = z_pad+z-8;

    x_end = x_start+8;
    y_end = y_start+8;
    z_end = z_start+8;
else
    x_start = x_pad+x-4;
    y_start = y_pad+y-4;
    z_start = z_pad+z;

    x_end = x_start+8;
    y_end = y_start+8;
    z_end = z_start+8;
end

mask = zeros(size(template));
mask(mm) = 1;

cor = zeros(size(target_image));
template1 = template(mm);
for zz = z_start:z_end
    for yy = y_start:y_end
        for xx = x_start:x_end
            tempx_start = xx-x_pad;
            tempx_end = xx+x_pad;
            tempy_start = yy-y_pad;
            tempy_end = yy+y_pad;
            tempz_start = zz-z_pad;
            tempz_end = zz+z_pad;
            
            temp = target_pad(tempx_start:tempx_end,tempy_start:tempy_end,...,
                              tempz_start:tempz_end);
            temp1 = temp(mm);
            
            corx = xx-x_pad;
            cory = yy-y_pad;
            corz = zz-z_pad;
            
            cor(corx,cory,corz) = cor_own(temp1,template1);
        end
    end
end