close all,clear all,clc;
%%
% 假设所有的.mat文件都存放在同一个文件夹中，并且每个文件中都只有一个列向量
folderPath = '/Users/taylor/Desktop/code/noise_stability_measurement/row'; % 存放.mat文件的文件夹路径%此向量记载了单幅496*496pixels全息图所有像素点的像素值
files = dir(fullfile(folderPath, '*.mat')); % 获取所有.mat文件

% 初始化一个空矩阵用于存储合并后的数据
% 假设你知道所有列向量的长度是相同的
numRows = size(load(fullfile(folderPath, files(1).name)), 2); % 获取列向量的长度
mergedData = []; % 初始化合并后的数据矩阵

% 循环读取每个文件并合并数据
for i = 1:length(files)
    matData = load(fullfile(folderPath, files(i).name)); % 加载.mat文件
    % 假设每个.mat文件中的变量名是相同的，例如'data'
    varName = fieldnames(matData); % 获取变量名
    mergedData = [mergedData, matData.(varName{1})]; % 合并列向量
end
% 将合并后的数据写入Excel文件
% outputExcelPath = '合并后的数据.xlsx'; % 输出Excel文件的路径和文件名
% xlswrite(outputExcelPath, mergedData); % 写入Excel文件

%% 
% 假设mergedData是你已经合并好的数据矩阵
% 确保mergedData至少有两列
if size(mergedData, 2) > 1
    % 将第一列作为基准值，从第二列到最后一列的每一列都减去第一列的对应元素
    newMatrix = mergedData(:, 2:end) - mergedData(:, 1);
else
    error('mergedData矩阵的列数少于两列，无法执行操作。');
end

%% Temporal noise sensitivity
rowStd = std(mergedData, 0, 2);
% 输出每行的标准差
disp('每行的标准差分别为：');
disp(rowStd);

% 假设newMatrix是你已经计算好的新矩阵
% 计算每行的均值
rowMeans = mean(newMatrix, 2);
% 计算每行的标准差
rowStds = std(newMatrix, 0, 2);
% 将每行的均值和标准差写入Excel文件
outputExcelPath = '行标准差.xlsx'; % 输出Excel文件的路径和文件名
% 写入每行的标准差
xlswrite(outputExcelPath, num2cell(rowStds), 'Sheet1', 'A1') % 从A2开始写入每行的标准差

%% Temporal noise sensitivity plot
% 绘制直方图
figure; % 创建一个新的图形窗口
histogram(rowStds, 'Normalization', 'count'); % 绘制直方图，并将直方图归一化为计数

% 添加标题和轴标签
%title('Distribution of Row Standard Deviations');
xlabel('Standard deviation(rad)');
ylabel('Number of counts');

% 可选：添加网格
grid off; % 添加网格线

% 可选：设置x轴和y轴的限制
% xlim([0 1]); % 如果需要限制x轴的范围，取消注释并替换minStd和maxStd为实际值
% ylim([minFreq maxFreq]); % 如果需要限制y轴的范围，取消注释并替换minFreq和maxFreq为实际值

%% Spartial noise sensitivity
% 假设newMatrix是你已经计算好的新矩阵
% 计算每列的标准差
colStds = std(newMatrix, 0, 1);

% 输出每列的标准差
disp('每列的标准差分别为：');
disp(colStds);

% 将每列的标准差写入Excel文件
outputExcelPath = '列标准差.xlsx'; % 输出Excel文件的路径和文件名

% 写入每列的标准差
xlswrite(outputExcelPath, num2cell(colStds), 'Sheet1', 'A1'); % 从A1开始写入每列的标准差

% 如果需要，可以添加更多的Excel文件操作，如设置列宽、添加标题等

%% Spartial noise sensitivity plot
% 绘制直方图
figure; % 创建一个新的图形窗口
histogram(colStds, 'Normalization', 'count'); % 绘制直方图，并将直方图归一化为计数

% 添加标题和轴标签
%title('Distribution of Column Standard Deviations');
xlabel('Standard deviation(rad)');
ylabel('Number of counts');

% 可选：添加网格
grid off; % 添加网格线

% 可选：设置x轴和y轴的限制
% xlim([minStd maxStd]); % 如果需要限制x轴的范围，取消注释并替换minStd和maxStd为实际值
% ylim([minFreq maxFreq]); % 如果需要限制y轴的范围，取消注释并替换minFreq和maxFreq为实际值

%% Spartialtemporal noise sensitivity
% 假设mergedData是你已经合并好的数据矩阵
% 将mergedData中的数据展平为一个列向量
dataVector = newMatrix(:);
% 计算所有元素的标准差
overallStd = std(dataVector);
% 输出标准差
disp(['整体标准差为：', num2str(overallStd)]);

% 创建一个表格，包含所有元素和标准差
output_table = table(dataVector, repmat(overallStd, length(dataVector), 1), ...
                     'VariableNames', {'DataVector', 'OverallStd'});

% 将表格写入新的Excel文件
output_filename = 'Data_and_Std.xlsx';
writetable(output_table, output_filename);
disp(['所有元素和整体标准差已输出到Excel文件：', output_filename]);

%% Spartialtemporal noise sensitivity plot 
% 假设newMatrix是你已经计算好的新矩阵
% 将newMatrix中的所有元素展平成一个列向量
allElements = newMatrix(:);

% 绘制直方图
figure; % 创建一个新的图形窗口
histogram(allElements, 'Normalization', 'count'); % 绘制直方图，并将直方图归一化为计数

% 添加标题和轴标签
%title('Distribution of All Elements in newMatrix');
xlabel('Value(rad)');
ylabel('Number of counts');

% 可选：添加网格
grid off; % 添加网格线

% 可选：设置x轴和y轴的限制
 xlim([-1 1]); % 如果需要限制x轴的范围，取消注释并替换minValue和maxValue为实际值
% ylim([minFreq maxFreq]); % 如果需要限制y轴的范围，取消注释并替换minFreq和maxFreq为实际值

% 可选：为直方图添加更多的定制，如颜色、条形宽度等
% histogram(allElements, 'Normalization', 'count', 'FaceColor', 'blue', 'BarWidth', 1);