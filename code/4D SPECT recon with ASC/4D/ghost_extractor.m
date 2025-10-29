% Interpolation 3D of the ghost information

function [Vx,Vy,Vz,heart_ref,heart_tar] = ghost_extractor(lgrX,lgrY,lgrZ,filename)

% INPUTS:
% [lgrX,lgrY,lgrZ]   : dimensions of the reference image corresponding of the
%                      motion file
% filename           : name of the motion file

% OUTPUTS:
% [Vx,Vy,Vz]         : heart motion matrices (same size of the reference
%                      image; motion of pixel not inside heart are set to 0) 
% heart_ref          : matrix with same size of the reference image; pixel 
%                      inside heart of the reference image is set to 1,
%                      otherwise 0
% heart_tar          : matrix with same size of the target image; pixel 
%                      inside heart of the target image is set to 1,
%                      otherwise 0
% ghost_database.mat : saving file of the outputs

%Extraction of the theorical values:
fid = fopen(filename);
fgetl(fid);
fgetl(fid);
C = textscan(fid,'%*s %*s %f %f %f %*s %f %f %f %*s %f %f %f');
fclose(fid);
Xref = C{1};
Yref = C{2};
Zref = C{3};
Xtar = C{4};
Ytar = C{5};
Ztar = C{6};
Vtx = C{7};
Vty = C{8};
Vtz = C{9};
[lgr,ind]=size(Xref);

% Interpolation are done under a circle of radius ctr_d
% Criterion of proximity
crt_d = 1;

%Initial estimates of the velocities are set as zero
Vx = zeros(lgrX,lgrY,lgrZ);
Vy = zeros(lgrX,lgrY,lgrZ);
Vz = zeros(lgrX,lgrY,lgrZ);
heart_ref = zeros(lgrX,lgrY,lgrZ);
heart_tar = zeros(lgrX,lgrY,lgrZ);

Xref = Xref + 1;
Yref = Yref + 1;
Zref = Zref + 1;
%Processing of the reference image
for i=1:lgrX
    for j=1:lgrY
        for k=1:lgrZ
            dtot=0;
            for l=1:lgr
                dx = abs(i-Xref(l));
                dy = abs(j-Yref(l));
                dz = abs(k-Zref(l));
                d = (dx^2+dy^2+dz^2)^(1/2);
                if d==0
                    heart_ref(i,j,k) = 1;
                    Vx(i,j,k) = Vtx(l);
                    Vy(i,j,k) = Vty(l);
                    Vz(i,j,k) = Vtz(l);
                    break;
                elseif d<crt_d
                    heart_ref(i,j,k) = 1;
                    dtot = dtot + 1/d;
                    Vx(i,j,k) = Vx(i,j,k) + (1/d)*Vtx(l);
                    Vy(i,j,k) = Vy(i,j,k) + (1/d)*Vty(l);
                    Vz(i,j,k) = Vz(i,j,k) + (1/d)*Vtz(l);
                end
            end
            if dtot~=0
                Vx(i,j,k) = Vx(i,j,k)/dtot;
                Vy(i,j,k) = Vy(i,j,k)/dtot;
                Vz(i,j,k) = Vz(i,j,k)/dtot;
            end
        end
    end
end

Xtar = Xtar + 1;
Ytar = Ytar + 1;
Ztar = Ztar + 1;
%Processing of the target image
for i=1:lgrX
    for j=1:lgrY
        for k=1:lgrZ
            for l=1:lgr
                dx = abs(i-Xtar(l));
                dy = abs(j-Ytar(l));
                dz = abs(k-Ztar(l));
                d = (dx^2+dy^2+dz^2)^(1/2);
                if d<crt_d
                    heart_tar(i,j,k) = 1;
                end
            end
        end
    end
end

%save('ghost_database','Vx','Vy','Vz','heart_ref','heart_tar');