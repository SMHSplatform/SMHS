% 假设你的Excel文件名为 'tSNE.xlsx'，并且数据从第一行第一列开始
filename = 'tSNE.xlsx';

% 读取Excel文件中的数据
data = readtable(filename);

% 提取特征数据（前81列）
features = table2array(data(:, 1:81));

% 提取细胞标签（第82列）
labels = data{:, 82};

% 进行t-SNE降维
Y = tsne(features, 'Algorithm', 'barneshut', 'Perplexity', 40, 'NumDimensions', 2);

% 为每种细胞类型指定一个颜色
colors = [
    0.49 0.18 0.56; % RBC - 紫色
    0.30 0.75 0.93; % PBMC - 蓝色
    0.47 0.67 0.19; % MDAMB231 - 绿色
];

% 绘制t-SNE降维后的聚类图
figure;
hold on; % 保持当前图像，允许在同一图像上绘制多个图层

% 绘制每种细胞类型的散点图
scatter(Y(strcmp(labels, 'RBC'), 1), Y(strcmp(labels, 'RBC'), 2), 10, colors(1,:), 'filled');
scatter(Y(strcmp(labels, 'PBMC'), 1), Y(strcmp(labels, 'PBMC'), 2), 10, colors(2,:), 'filled');
scatter(Y(strcmp(labels, 'MDAMB231'), 1), Y(strcmp(labels, 'MDAMB231'), 2), 10, colors(3,:), 'filled');

hold off; % 释放图像

title('t-SNE Clustering');
xlabel('t-SNE Feature 1');
ylabel('t-SNE Feature 2');

% 创建图例标签
legend_labels = {'RBC', 'PBMC', 'MDA-MB-231'};
legend(legend_labels, 'Location', 'best');