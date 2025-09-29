% Program to analyse and visualize distribution of percent of 1s in SRAM power-up states for Self-Referencing Tests for Counterfeit Detection.

%% =========================================================
%  SCRIPT: Multi-Chip Power-Up State Analysis
%  Author: Zakia Tisha
%  Email: zakia.tisha@auburn.edu
%  Institute: Auburn University
%  Advisor: Dr. Ujjwal Guin
%  PURPOSE:
%   - Read CSV files for single chips or chip pairs
%   - Odd chip → 1…N states, Even chip → N+1…2N states
%   - Filter data (%1s within thresholds)
%   - Fit Normal distributions, extract sigma values
%   - Plot sigma vs days for user-chosen block sizes
%   - Plot sigma vs block sizes (Day 1 vs Final Day)
% ===========================================================

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

Days = (0:n-1);

%% ============================================================
% STEP 3: Fit Normal Distributions
% ============================================================
activeChips = find(~cellfun(@isempty,ChipData));
FitResults = cell(max(activeChips), n);

for chip = activeChips
    for day = 1:n
        FitResults{chip,day} = cell(size(blockSizes));
        for idx = 1:length(blockSizes)
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
% STEP 4: Extract Sigma Values
% ============================================================
SigmaValues = cell(max(activeChips),1);

for chip = activeChips
    SigmaMatrix = nan(length(blockSizes), n);
    for day = 1:n
        for idx = 1:length(blockSizes)
            if ~isempty(FitResults{chip,day}{idx})
                SigmaMatrix(idx,day) = FitResults{chip,day}{idx}.sigma;
            end
        end
    end
    SigmaValues{chip} = SigmaMatrix;
end

%% ============================================================
% STEP 5: Plot Results
% ============================================================

% --- User selects block sizes for Sigma vs Days ---
disp('Available block sizes:');
disp(blockSizes);
selectedBlocks = input('Enter block sizes as a vector (e.g., [64 128]): ');
if any(~ismember(selectedBlocks, blockSizes))
    error('One or more selected block sizes are not in the available list!');
end

% === Sigma vs Days (combined plot for all selected block sizes) ===
figure;
hold on;

markers = {'-o','-s','-d','-^','-v','-p','-h'}; % marker styles
colors  = lines(max(activeChips));             % color set for chips

legendEntries = {};
plotHandles = [];

for b = 1:length(selectedBlocks)
    block = selectedBlocks(b);
    idx   = find(blockSizes == block);

    for chip = activeChips
        if isempty(ChipData{chip}), continue; end

        % Build the sigma vector across days
        sigmaVals = nan(1,n);
        for day = 1:n
            fitStruct = FitResults{chip,day}{idx};
            if ~isempty(fitStruct)
                sigmaVals(day) = fitStruct.sigma;
            end
        end

        % Plot
        markerStyle = markers{mod(b-1, numel(markers)) + 1};
        h = plot(Days, sigmaVals, markerStyle, ...
            'LineWidth',1.5, 'Color', colors(chip,:));

        % Add legend entry only if data exists
        if any(~isnan(sigmaVals))
            plotHandles(end+1) = h;
            legendEntries{end+1} = sprintf('Chip %d (%dx%d)', ...
                chip, selectedBlocks(b), selectedBlocks(b));
        end
    end
end

grid on; box on;
set(gca,'FontSize',14)
xlabel('Days','FontSize',16)
ylabel('\sigma','FontSize',20)
legend(plotHandles, legendEntries, 'FontSize', 12, 'Location','best');
title('Sigma vs Days (all block sizes)','FontSize',16);
hold off;

% --- Sigma vs Block Sizes (Day 1 vs Final Day) ---
figure;
hold on;
for chip = activeChips
    plot(blockSizes, SigmaValues{chip}(:,1), '-o', 'LineWidth',1.5);
    plot(blockSizes, SigmaValues{chip}(:,end), '-d', 'LineWidth',1.5);
end
grid on; box on;
set(gca,'FontSize',16);
xticks([12 32 42 51 64 85 128]);
xticklabels({'12x12','32x32','42x42','51x51','64x64','85x85','128x128'});
xlabel('Block Sizes','FontSize',16);
ylabel('\sigma','FontSize',23);
legendStr = [arrayfun(@(x) sprintf('Day1 C%d',x),activeChips,'UniformOutput',false), ...
             arrayfun(@(x) sprintf('Day%d C%d',n,x),activeChips,'UniformOutput',false)];
legend(legendStr,'FontSize',12);
title('Sigma vs Block Sizes (Day 1 vs Final Day)');






