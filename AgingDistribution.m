clc;
clear all

% Prompt the user to enter the number of datasets (n)
fileindex1 = input('enter the number of power up states for one Chip: ');
fileindex2 = 100;
fileindex3 = fileindex1;

n=2;

% Define the range for the second column
minRange = 0.55; % Replace with your desired minimum value
maxRange = 0.74 ; % Replace with your desired maximum value

%%Loop through and read the specified CSV files for each set of data
for i = 1:n
    %%Prompt the user to enter the base filename for each dataset
    baseFilename1 = input(['Enter the base filename for Chip 1 and Chip 2- dataset ' num2str(i) ': '], 's');

    for j = 1:fileindex1
        %%Construct the file name for the current data set
        filename1 = sprintf('%s%dAvgList.csv', baseFilename1, j + fileindex3);

        %%Check if the file exists
        if exist(filename1, 'file') == 2
            %%Read the data if the file exists
            data1 = readmatrix(filename1);
            
            %%Filter the data based on the range for the second column
            validIndices1 = data1(:, 2) >= minRange & data1(:, 2) <= maxRange;
            filteredData1 = data1(validIndices1, :);
            filteredData1(:, 2) = filteredData1(:, 2) * 100;
            
            %%Store the filtered data
            Chip1_read{j} = filteredData1;

             else
            fprintf('File not found: %s\n', filename1);
        end
            
        %%Construct the file name for the current data set
        filename2 = sprintf('%s%dAvgList.csv', baseFilename1, j + fileindex3);

        %%Check if the file exists
        if exist(filename2, 'file') == 2
            %%Read the data if the file exists
            data2 = readmatrix(filename2);
            
            %%Filter the data based on the range for the second column
            validIndices2 = data2(:, 2) >= minRange & data2(:, 2) <= maxRange;
            filteredData2 = data2(validIndices2, :);
             filteredData2(:, 2) = filteredData2(:, 2) * 100;
            
            %%Store the filtered data
            Chip2_read{j} = filteredData2;
        else
            fprintf('File not found: %s\n', filename);
        end
    end

    %%Store the data for this dataset
    Chip1Data{i} = Chip1_read;
    Chip2Data{i} = Chip2_read;
end



%Initializing block sizes and conditions
conditions = {4, 9, 16, 25, 36, 49, 64, 81,100, 121, 144, 169, 196, 225, 256, 289, 324, 361, 400};
sizes = {128, 85, 64, 51, 42, 36, 32, 28, 25, 23, 21, 19, 18, 17, 16, 15, 14, 13, 12};


blockSizes = [128, 85, 64, 51, 42, 36, 32, 28, 25, 23, 21, 19, 18, 17, 16, 15, 14, 13, 12];

%blockSizes = [85, 64, 51, 42, 32];


for chipIdx = 1:2
    %Dynamically construct variable names
    chipDataVarName = sprintf('Chip%dData', chipIdx);
    chipSigmaVarName = sprintf('Chip%d_Sigma', chipIdx);
    chipTVarName = sprintf('Chip%d_t', chipIdx);

    %Access the data using eval
    ChipData = eval(chipDataVarName);
    
    for i = 1:n
        for j = 1:length(conditions)
            for k = 1:numel(ChipData{i})
                condition = (ChipData{i}{k}(:, 1) == conditions{j});
                eval(sprintf('%s{%d, %d, %d} = ChipData{i}{k}(condition, 2);', chipSigmaVarName, j, i, k));
            end

            %%Concatenating the data vertically
            eval(sprintf('%s{%d, %d} = vertcat(%s{%d, %d, :});', chipTVarName, j, i, chipSigmaVarName, j, i));
        end
    end
end

%Access the final results
Chip1 = Chip1_t.';
Chip2 = Chip2_t.';



blockSizes = [128, 85, 64, 51, 42, 36, 32, 28, 25, 23, 21, 19, 18, 17, 16, 15, 14, 13, 12];

Days = (0:n-1);

for Z = blockSizes
    for i = 1:n
         eval(['Chip1_Day' num2str(i) '_' num2str(Z) 'x' num2str(Z) ' = cell2mat(Chip1(i,' num2str(find(blockSizes==Z)) '));']);
         eval(['Chip2_Day' num2str(i) '_' num2str(Z) 'x' num2str(Z) ' = cell2mat(Chip2(i,' num2str(find(blockSizes==Z)) '));']);
    end
end



numBins = 20;



for i = 1:n
    for Z = blockSizes
        eval(['Chip1_Day' num2str(i) '_fitdists = cell(size(blockSizes));']);
        eval(['Chip2_Day' num2str(i) '_fitdists = cell(size(blockSizes));']);
    end
end

for i = 1:n
    for Z = blockSizes
        eval(['Chip3_Day' num2str(i) '_fitdists = cell(size(blockSizes));']);
        eval(['Chip4_Day' num2str(i) '_fitdists = cell(size(blockSizes));']);
    end
end

for idx = 1:length(blockSizes)
    X = blockSizes(idx);

    for day = 1:n  % Loop through Day1, Day2, ..., DayN
        varName1 = ['Chip1_Day' num2str(day) '_fitdist'];
        Chip1_fitDist = fitdist(eval(['Chip1_Day' num2str(day) '_' num2str(X) 'x' num2str(X)]), 'Normal');
        assignin('base', varName1, Chip1_fitDist);
        % %Store fitdist for each day in separate variables
       eval(['Chip1_Day' num2str(day) '_fitdists{' num2str(idx) '} = Chip1_fitDist;']);


        varName2 = ['Chip2_Day' num2str(day) '_fitdist'];
        Chip2_fitDist = fitdist(eval(['Chip2_Day' num2str(day) '_' num2str(X) 'x' num2str(X)]), 'Normal');
        assignin('base', varName2, Chip2_fitDist);
        % %Store fitdist for each day in separate variables
        eval(['Chip2_Day' num2str(day) '_fitdists{' num2str(idx) '} = Chip2_fitDist;']);

    end
end

% Assuming you have a loop index 'i'
n_values = zeros(2, 1); % Initialize an array to store n1, n2, n3

for i = 1:2
      
    day_data = eval(['Chip1_Day' num2str(i) '_64x64']); % Evaluate the variable name dynamically
    n_values(i) = (max(day_data(:)) - min(day_data(:))) / 0.008;
end

n1 = ceil(n_values(1)/100);
n2 = ceil(n_values(2)/100);

figure
histfit(Chip1_Day1_64x64, n1)
hold on
histfit(Chip1_Day2_64x64, n2)

set(gca,'FontSize',10)
xlabel('p1s','FontSize',12);
ylabel('Frequency','FontSize', 12);
%title('Block Size 64x64');
%legend('New Chip',' ','7 days aged chip', ' ','FontSize', 20, 'fontname','Times New Roman')

