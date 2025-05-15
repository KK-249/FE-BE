clear
clc

status = system('FinalBEM.exe');

%—— 1. 确认当前目录和文件存在 ——————————————————
cwd = pwd;
fprintf('MATLAB 工作目录：%s\n', cwd);
if ~exist('FinalBEM.exe','file')
    error('当前目录里找不到 FinalBEM.exe！');
end
if ~exist('example66.dat','file')
    error('当前目录里找不到 example66.dat！');
end

%—— 2. 用 cmd /C 先切目录再运行，并把所有输出重定向到日志 ——
logfile = fullfile(cwd, 'run.log');
cmd = sprintf('cmd /C "cd /D %s && FinalBEM.exe > "%s" 2>&1"', cwd, logfile);
status = system(cmd);

%—— 3. 检查返回码和日志 ——————————————————————
if status ~= 0
    fprintf('>> FinalBEM.exe 返回错误码 %d ，请打开 %s 查看详细信息：\n', status, logfile);
    type(logfile)
else
    fprintf('>> FinalBEM.exe 运行完毕（错误码 %d），日志在 %s：\n', status, logfile);
    type(logfile)
end

%—— 4. 验证输出文件是否生成 ————————————————————
if exist('output.dat','file')
    fprintf('>>> output.dat 已生成 (%d 字节)\n', dir('output.dat').bytes);
else
    warning('>>> output.dat 未生成，请检查 run.log 中的报错！');
end