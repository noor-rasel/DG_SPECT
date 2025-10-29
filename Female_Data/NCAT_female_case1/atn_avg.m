%% === Average attenuation maps across 8 respiratory phases ===
path_setting1;   % if you already define atn_act_data_root_folder
atn_data = atn_act_data_root_folder;

num_resp = 8;             
nx = 64; ny = 64; nz = 64;

% Initialize array to store summed attenuation maps
atn_sum = zeros(nx, ny, nz, 'single');

for resp = 1:num_resp
    % Build full file path
    filename = sprintf('%s/cardiac8_atn_%d.bin', atn_data, resp);
    
    % Open file
    fid = fopen(filename, 'rb');
    if fid == -1
        error('Cannot open file: %s', filename);
    end
    
    % Read data
    tmp = fread(fid, nx*ny*nz, 'single');
    fclose(fid);

    % Reshape to 3D matrix
    tmp = reshape(tmp, [nx, ny, nz]);

    % Apply same orientation as other attenuation maps
    tmp = flip(tmp, 1);
    tmp = flip(tmp, 3);
    tmp = rot90(tmp, -1);
    
    % Accumulate
    atn_sum = atn_sum + tmp;
    fprintf('Loaded: %s\n', filename);
end

% Compute average attenuation map
attn_avg = atn_sum / num_resp;

% Save averaged attenuation map
output_filename = fullfile(atn_data, 'cardiac8_atn_avg.bin');
fid = fopen(output_filename, 'wb');
fwrite(fid, attn_avg, 'single');
fclose(fid);

fprintf('Saved averaged attenuation map → %s\n', output_filename);

% === Optional visualization ===
figure;
imagesc(attn_avg(:,:,41));
axis image off;
colormap(gray);
title('Averaged Attenuation Map (Slice 41)');

