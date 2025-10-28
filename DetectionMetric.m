% ============================================================
% Script to Compute Δ = σ_51x51 − σ_128x128 for new and aged chips
%
% Author: Zakia Tisha
% Advisor: Dr. Ujjwal Guin
% Institute: Auburn University
% ============================================================

clc; clear all;

%% -----------------------------
% User Inputs
% -----------------------------
chipChoice = input('Do you want to analyze (1) Single chip(s) or (2) Chip pairs? Enter 1 or 2: ');
fileindex1 = input('Enter the number of power-up states per chip (e.g., 100): ');
n          = input('Enter the number of datasets (days) (e.g., 8): ');

% Filtering thresholds
minRange = 0.5;
maxRange = 0.85;

% Containers
ChipData = cell(1, 8); % up to 8 chips supported

%% ============================================================
% STEP 1: Read and Filter CSV Data
% ============================================================
if chipChoice == 1
    %% --- Single Chip Mode ---
    selectedChips = input('Enter chip number(s) to analyze individually (e.g., [1], [2], [3 5]): ');

    for chip = selectedChips
        for day = 1:n
            datasetDir = input(['Enter dataset directory for Chip ' num2str(chip) ...
                                ' - Dataset ' num2str(day) ': '], 's');

            % Ensure underscore at the end
            if datasetDir(end) ~= '_'
                datasetDir = [datasetDir '_'];
            end

            chip_read = cell(1,fileindex1);

            if mod(chip,2)==1
                % Odd chip → states 1…N
                for j = 1:fileindex1
                    filename = sprintf('%s%dAvgList.csv', datasetDir, j);
                    if exist(filename,'file')==2
                        data = readmatrix(filename);
                        valid = data(:,2) >= minRange & data(:,2) <= maxRange;
                        chip_read{j} = data(valid,:);
                    else
                        fprintf('File not found: %s\n', filename);
                        chip_read{j} = [];
                    end
                end
            else
                % Even chip → states N+1…2N
                for j = 1:fileindex1
                    filename = sprintf('%s%dAvgList.csv', datasetDir, j+fileindex1);
                    if exist(filename,'file')==2
                        data = readmatrix(filename);
                        valid = data(:,2) >= minRange & data(:,2) <= maxRange;
                        chip_read{j} = data(valid,:);
                    else
                        fprintf('File not found: %s\n', filename);
                        chip_read{j} = [];
                    end
                end
            end
            ChipData{chip}{day} = chip_read;
        end
    end

elseif chipChoice == 2
    %% --- Chip Pair Mode ---
    numPairs = input('Enter the number of chip pairs to include (e.g., 3 for Chips 1–6): ');

    for pair = 1:numPairs
        oddChip  = 2*pair - 1;
        evenChip = 2*pair;

        for day = 1:n
            datasetDir = input(['Enter dataset directory for Chip ' num2str(oddChip) ...
                                ' & Chip ' num2str(evenChip) ...
                                ' - Dataset ' num2str(day) ': '], 's');

            % Ensure underscore at the end
            if datasetDir(end) ~= '_'
                datasetDir = [datasetDir '_'];
            end

            % Odd chip (1…N)
            chip_read = cell(1,fileindex1);
            for j = 1:fileindex1
                filename = sprintf('%s%dAvgList.csv', datasetDir, j);
                if exist(filename,'file')==2
                    data = readmatrix(filename);
                    valid = data(:,2) >= minRange & data(:,2) <= maxRange;
                    chip_read{j} = data(valid,:);
                else
                    fprintf('File not found: %s\n', filename);
                    chip_read{j} = [];
                end
            end
            ChipData{oddChip}{day} = chip_read;

            % Even chip (N+1…2N)
            chip_read = cell(1,fileindex1);
            for j = 1:fileindex1
                filename = sprintf('%s%dAvgList.csv', datasetDir, j+fileindex1);
                if exist(filename,'file')==2
                    data = readmatrix(filename);
                    valid = data(:,2) >= minRange & data(:,2) <= maxRange;
                    chip_read{j} = data(valid,:);
                else
                    fprintf('File not found: %s\n', filename);
                    chip_read{j} = [];
                end
            end
            ChipData{evenChip}{day} = chip_read;
        end
    end
end

%% ============================================================
% STEP 2: Setup Conditions and Block Sizes
% ============================================================
conditions = {4, 9, 16, 25, 36, 49, 64, 81, 100, 121, 144, ...
              169, 196, 225, 256, 289, 324, 361, 400};
blockSizes = [128, 85, 64, 51, 42, 36, 32, 28, 25, ...
              23, 21, 19, 18, 17, 16, 15, 14, 13, 12];

idx_128 = find(blockSizes == 128);
idx_51  = find(blockSizes == 51);

%% ============================================================
% STEP 3: Fit Normal Distributions
% ============================================================
activeChips = find(~cellfun(@isempty,ChipData));
FitResults = cell(max(activeChips), n);

for chip = activeChips
    for day = 1:n
        FitResults{chip,day} = cell(size(blockSizes));
        for idx = [idx_51, idx_128]  % only process 51×51 and 128×128
            data_cat = [];
            for k = 1:numel(ChipData{chip}{day})
                data_k = ChipData{chip}{day}{k};
                if ~isempty(data_k)
                    mask = (data_k(:,1) == conditions{idx});
                    data_cat = [data_cat; data_k(mask,2)];
                end
            end
            if ~isempty(data_cat)
                FitResults{chip,day}{idx} = fitdist(data_cat,'Normal');
            else
                FitResults{chip,day}{idx} = [];
            end
        end
    end
end

%% ============================================================
% STEP 4: Compute Δ for each chip (New vs Aged)
% ============================================================
fprintf('\n=== Δ (Sigma_51x51 − Sigma_128x128) for New and Aged Chips ===\n\n');
fprintf('%-8s %-15s %-15s\n','Chip','# New (Day 0)','Aged (Final Day)');

for chip = activeChips
    sigma_new_51  = FitResults{chip,1}{idx_51}.sigma;
    sigma_new_128 = FitResults{chip,1}{idx_128}.sigma;
    sigma_old_51  = FitResults{chip,n}{idx_51}.sigma;
    sigma_old_128 = FitResults{chip,n}{idx_128}.sigma;

    delta_new = sigma_new_51 - sigma_new_128;
    delta_old = sigma_old_51 - sigma_old_128;

    fprintf('Chip %-3d  %-15.3f %-15.3f\n', chip, delta_new, delta_old);
end

fprintf('\nComputation complete. Δ values printed above.\n');
