clear
clc

% === 参数设置 ===
   
dataFolder = 'C:\Users\Administrator\Desktop\有限元程序\result';
% dataFile   = 'C:\Users\Administrator\Desktop\有限元程序\result\output_ 60.dat';   % 要读取的文件
lineRange  = 20:25;            % 要读取的行号
colIndex   = 11;               % 要读取的列号

A1 = zeros(1000,1);A2 = zeros(1000,1);A3 = zeros(1000,1);
A4 = zeros(1000,1);A5 = zeros(1000,1);A6 = zeros(1000,1);

for k = 1:1000
    if k <= 99 
    fname  = sprintf('output_  %2d.dat', k);
    dataFile = fullfile(dataFolder, fname);
    % === 读入整个文件为一个大字符向量 ===
    txtAll = fileread(dataFile);
    elseif k ==1000
    fname  = sprintf('output_%2d.dat', k);
    dataFile = fullfile(dataFolder, fname);
    % === 读入整个文件为一个大字符向量 ===
    txtAll = fileread(dataFile);   
    else
     fname  = sprintf('output_ %2d.dat', k);
    dataFile = fullfile(dataFolder, fname);
    % === 读入整个文件为一个大字符向量 ===
    txtAll = fileread(dataFile);   
    end

% === 按行拆分 ===
lines = regexp(txtAll, '\r\n|\r|\n', 'split');

% === 检查是否足够行数 ===
if numel(lines) < max(lineRange)
    error('文件行数不足 %d 行！', max(lineRange));
end

% === 预分配存放第17列的值 ===
nLines = numel(lineRange);
col11_vals = nan(nLines,1);

% === 遍历 12 到 17 行，提取第 17 列 ===
for i = 1:nLines
    ln = strtrim(lines{ lineRange(i) });
    % 按空格、逗号、制表符分隔字段
    tokens = regexp(ln, '[\s,]+', 'split');

    % 2) 去掉可能的空字符串
    tokens(cellfun('isempty', tokens)) = [];
%     tokens = regexp(ln, '[, \t]+', 'split');
    % 容错：检查列数
    if numel(tokens) < colIndex
        fprintf('第 %d 行被拆成 %d 个字段\n', i, numel(tokens));
        error('第 %d 行只有 %d 列，无法取到第 %d 列', lineRange(i), numel(tokens), colIndex);
    end
    % 转为数字
    val = str2double(tokens{colIndex});
    if isnan(val)
        error('第 %d 行第 %d 列不是一个有效数字："%s"', ...
              lineRange(i), colIndex, tokens{colIndex});
    end
    col11_vals(i) = val;
end

% 如果你真的需要把它们存成 6 个单独的数组（或变量），可以这样：
A1(k,1) = col11_vals(1);
A2(k,1) = col11_vals(2);
A3(k,1) = col11_vals(3);
A4(k,1) = col11_vals(4);
A5(k,1) = col11_vals(5);
A6(k,1) = col11_vals(6);
end

A1(:,1) = 20*log10(A1(:,1)./sqrt(2)*1000000);    %声压级总是 总声压的均方根值
A2(:,1) = 20*log10(A2(:,1)./sqrt(2)*1000000);
A3(:,1) = 20*log10(A3(:,1)./sqrt(2)*1000000);
A4(:,1) = 20*log10(A4(:,1)./sqrt(2)*1000000);
A5(:,1) = 20*log10(A5(:,1)./sqrt(2)*1000000);
A6(:,1) = 20*log10(A6(:,1)./sqrt(2)*1000000);

plot(1:1000,A1(1:1000,1));
plot(1:1000,A2(1:1000,1));
plot(1:1000,A3(1:1000,1));
plot(1:1000,A4(1:1000,1));
plot(1:1000,A5(1:1000,1));
plot(1:1000,A6(1:1000,1));