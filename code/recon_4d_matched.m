path_setting;

addpath(recon_root_folder+"code/");       
addpath(recon_root_folder+"code/4D SPECT recon with ASC/4D/");
addpath(recon_root_folder+"code/4D SPECT recon with ASC/weight_gen/");
addpath(recon_root_folder+"code/extra_files/",'-begin');

data_root_folder = simind_attn_data_root_folder + "all_phases";
projections = zeros(64,64,64,8);
recon_3d_card_phases = zeros(64,64,64,8,8);
recon_4d_card_phases = zeros(64,64,64,8,8);
recon_4d_card_resp_phases = zeros(64,64,64,8,8);
total_counts = 8e6;
total_counts_rp_1 = total_counts / 8;
blur = 1; 
OF_tag=0; 
load('roi.mat');

% === Get seed from environment ===
rand_seeds = [1, 4477, 8980, 1235321, 888888888, ...
             19971212, 123456789, 84847234, 6950321, ...
             123, 999999999, 7832645, 5543689, ...
             77234441, 2869, 13579, 595964237, ...
             2368475, 39645, 47, 10101, 8764101, ...
             4120339, 249988, 22740, 589033, ...
             657723, 20230110, 687321365, 95632113];

seed_index = str2double(getenv('seed_i'));
rng(rand_seeds(seed_index));  % Seed the random number generator

% rng(20230110);

for resp = 0 : 7
    % fprintf('\n>>> Respiratory Phase %d <<<\n', resp + 1);

    %% === Load Simind Projections Data for All Cardiac Phases ===
    for cardica = 0 : 7
        filename = sprintf('%s/cardica%d_resp_%d/cardiac%d_tot_w1.a00', data_root_folder, cardica+1, resp+1, cardica+1);
        fid = fopen(filename,'rb');
        simind = fread(fid,'float');
        fclose(fid);
        simind = reshape(simind,64,64,64);        
        projections(:,:,:,cardica+1) = simind;
    end

    %% === Rotate projections and normalize ===
    for i = 1:8
        cur_proj = projections(:,:,:,i);
        for j = 1:64
            cur_proj(:,:,j) = rot90(cur_proj(:,:,j),3);
        end
        projections(:,:,:,i) = cur_proj;
    end

    % fprintf('Total counts before norm: %d\n', sum(projections(:)));
    projections=projections/sum(projections(:)) * total_counts_rp_1;
    % fprintf('Total counts after norm: %d', sum(projections(:)));
    projections = random('poiss', projections);

    %% preparing system matrix using attenuation map
    load weight64_mn        
    filename = sprintf('%s/cardiac8_atn_%d.bin', simind_attn_data_root_folder, resp+1); % matched attenuation map
    fid=fopen(filename,'r');
    x=fread(fid,64^3,'single');
    fclose(fid);
    x=reshape(x,[64 64 64]);
    x = flip(x, 1);
    x = flip(x, 3); 
    wp_attnwgt=cell(4096,1); 

    % calculation of attenuation correction factors
    for n=1:64
        if n<33
            wp_attnwgt((n-1)*64+1:n*64)=attnmat(x,n,wp_vray((n-1)*64+1:n*64),...
                wp_ipxl((n-1)*64+1:n*64),wp_wgt((n-1)*64+1:n*64));  

        else 
            m=n-32;
            wp_attnwgt((n-1)*64+1:n*64)=attnmat(x,n,wp_vray((m-1)*64+1:m*64),...
                wp_ipxl((m-1)*64+1:m*64),wp_wgt((m-1)*64+1:m*64));
        end
    end

    % Find maximum attenuation factor per angle
    for j=1:64
        for i=1:64
            wp=wp_attnwgt{(j-1)*64+i}; 
            wp_S(i)=max(sum(wp)); 
        end
        wp_M(j)=max(wp_S); 
    end 

    % find maximum attenuation factor across all angles and normalize for all detectors 
    for j=1:64
        for i=1:64
            wp_attnwgt{(j-1)*64+i}=wp_attnwgt{(j-1)*64+i}/max(wp_M);
        end
    end 
    save('weight64_attn1.mat', 'wp_attnwgt');

    %% === 3D Reconstruction ===
    gbeta=0; 
    sbeta=7e-3; 
    sub_num=16; 
    it_num=10; 
    Im_maps=mbsrem4dv2(projections,repmat(roi,[1,1,64,8]),sub_num,it_num,OF_tag,...
         sbeta,gbeta,blur,1,0,0);
    % fprintf('3D noisy image data counts before norm: %d\n', sum(Im_maps(:)));
    Im_maps=Im_maps/sum(Im_maps(:)) * total_counts_rp_1;
    % fprintf('3D noisy image data counts after norm: %d\n', sum(Im_maps(:)));
    recon_3d_card_phases(:,:,:,:,resp+1) = Im_maps;
    save('recon_3d_card_phases_matched.mat', "recon_3d_card_phases");

    %% === Estimate Cardiac Motion===
    G = 8;
    yin=zeros(size(Im_maps));
    for g=1:G
        yin(:,:,:,g)=convn(Im_maps(:,:,:,g),ones(3,3,3)/27,'same');
    end

    M1=zeros(64,64,64,3,G);
    for i=1:G-1
        [vx,vy,vz]=motionele3d(yin(:,:,:,i:i+1),50,0.3);
        M1(:,:,:,:,i)=cat(4,vx,vy,vz);
    end
    [vx,vy,vz]=motionele3d(yin(:,:,:,[G 1]),50,0.3);
    M1(:,:,:,:,G)=cat(4,vx,vy,vz);

    for g=1:G
        ind5(g,:)=mod((g-2:g+2),G);
    end
    ind5(ind5==0)=G;

    for gate=1:G
        MM=mf2matrix3d_5g(M1(:,:,:,:,ind5(gate,:)),1);
        save(sprintf('n4dMM%d.mat', gate), 'MM');
    end

    %% === 4D Reconstruction ===
    gbeta=6e-3; 
    sbeta=0; 
    sub_num=16;
    it_num=10;
    Im_mapst=mbsrem4dv2(projections,repmat(roi,[1,1,64,8]),sub_num,it_num,OF_tag,...
        sbeta,gbeta,blur,1,0,0);
    % fprintf('4D noisy image data counts before norm: %d\n', sum(Im_mapst(:)));
    Im_mapst=Im_mapst/sum(Im_mapst(:)) * total_counts_rp_1;
    % fprintf('4D noisy image data counts after norm: %d\n', sum(Im_mapst(:)));
    recon_4d_card_phases(:,:,:,:,resp+1) = Im_mapst;
    save('recon_4d_card_phases_matched.mat', "recon_4d_card_phases");

    % === Clean up Temp Files ===
    fprintf('Cleaning up temporary files for respiratory phase %d\n', resp + 1);
    delete('weight64_attn1.mat');
    for gate = 1:8
        filename = sprintf('n4dMM%d.mat', gate);
        if exist(filename, 'file')
            delete(filename);
        end
    end
end 

%% === Estimate Respiratory Motion ===
recon_4d_card_phases_reduced = squeeze(sum(recon_4d_card_phases, 4)); % sum across card phases
recon_4d_card_phases_DVF = zeros(64,64,64,3,8,8);

for i = 1: 8
    fixed = recon_4d_card_phases_reduced(:,:,:,i);
    for j = 1: 8
        if(i==j), continue; end
        moving = recon_4d_card_phases_reduced(:,:,:,j);
        [DVF,movingReg] = imregdemons(moving,fixed,200,...
        'AccumulatedFieldSmoothing',1.3, 'DisplayWaitbar',false);
        recon_4d_card_phases_DVF(:,:,:,:,i,j) = DVF;
    end
end

%% === Temporal Smoothing using DVF ===
for i = 1:8
    for j = 1:8
        recon_4d_card_resp_phases(:,:,:,i,j) = recon_4d_card_phases(:,:,:,i,j);
        for k = 1 : 8
            if (k == j), continue; end               
            moving = recon_4d_card_phases(:,:,:,i,k);            
            DVF = recon_4d_card_phases_DVF(:,:,:,:,j,k);
            warped_image = imwarp(moving, DVF);
            weight = abs(1 - 2*abs(k-j)/8);
            recon_4d_card_resp_phases(:,:,:,i,j) = recon_4d_card_resp_phases(:,:,:,i,j) ...
                                             + weight * warped_image;
        end
    end
end

save('recon_4d_card_resp_phases_matched.mat',"recon_4d_card_resp_phases");

