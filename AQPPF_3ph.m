% AQPPF_3PH.M
% 基于随机变量二次型的三相不平衡配电网概率分布解析计算方法
% 参考：主要内容.md（第1-3章）+ 算例分析.md（第4章）
%
% 版本: v1.0
% 日期: 2026-02-07

dbstop if error; clc; clear; close all;

%% ========================================================================
%  段落1: 初始化
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  AQPPF_3ph - 三相不平衡配电网概率分布解析计算\n');
fprintf('%s\n\n', repmat('=',1,70));

% --- 路径配置 ---
projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'gncs'));
addpath(projectDir);  % 使 +tools 命名空间包可通过 tools.xxx 调用

% --- 负荷数据 ---
load(fullfile(projectDir, 'data', 'Load_Data.mat'));  % -> load_data (table)

% --- OpenDSS 算例配置 ---
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 修改这里
caseName = 'IEEE13';  % 可选: 'IEEE13', 'IEEE123'
gmm_num  = 20;        % GMM 分量数 K
sample_num = 10000;   % GMM 采样数
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch caseName
    case 'IEEE13'
        dssFile = 'C:\Program Files\OpenDSS\IEEETestCases\13Bus\IEEE13Nodeckt.dss';
    case 'IEEE123'
        dssFile = 'C:\Program Files\OpenDSS\IEEETestCases\123Bus\IEEE123Master.dss';
    otherwise
        error('不支持的算例: %s', caseName);
end
assert(isfile(dssFile), '[依赖缺失] OpenDSS算例文件不存在: %s', dssFile);

% --- 创建输出目录 ---
resultDir = fullfile(projectDir, 'result', caseName);
cacheDir  = fullfile(projectDir, 'cache');
if ~isfolder(resultDir), mkdir(resultDir); end
if ~isfolder(cacheDir),  mkdir(cacheDir);  end

fprintf('[配置] 算例: %s\n', caseName);
fprintf('[配置] GMM分量数 K=%d, 采样数=%d\n', gmm_num, sample_num);
fprintf('[配置] OpenDSS文件: %s\n', dssFile);

%% ========================================================================
%  段落2: 蒙特卡洛模拟（OpenDSS三相潮流，带缓存）
%  ========================================================================
N = height(load_data);  % 样本数量
cacheFile = fullfile(cacheDir, sprintf('MC_3ph_%s_T%d.mat', caseName, N));

if isfile(cacheFile)
    fprintf('\n[缓存] 发现缓存文件: %s\n', cacheFile);
    fprintf('[缓存] 加载缓存数据...\n');
    cacheData = load(cacheFile);
    V_mc_all     = cacheData.V_mc_all;
    nodeNames    = cacheData.nodeNames;
    nc           = cacheData.nc;
    Y3ph         = cacheData.Y3ph;
    nodeOrderY   = cacheData.nodeOrderY;
    loadNames    = cacheData.loadNames;
    lineNames    = cacheData.lineNames;
    linePower_mc = cacheData.linePower_mc;
    nLine        = cacheData.nLine;
    busNames     = cacheData.busNames;
    t1 = 0;
    fprintf('[缓存] 加载完成\n');

    % 重建负荷信息（对比方法需要）
    DSSObj_tmp = actxserver('OpenDSSEngine.DSS');
    DSSObj_tmp.Start(0);
    DSSObj_tmp.Text.Command = ['Compile "', dssFile, '"'];
    DSSCircuit_tmp = DSSObj_tmp.ActiveCircuit;
    iLoad_tmp = DSSCircuit_tmp.Loads;
    nLoad = iLoad_tmp.Count;
    loadKW_base  = zeros(nLoad, 1);
    loadKvar_base = zeros(nLoad, 1);
    idx_tmp = iLoad_tmp.First; li_tmp = 1;
    while idx_tmp > 0
        loadKW_base(li_tmp)  = iLoad_tmp.kW;
        loadKvar_base(li_tmp) = iLoad_tmp.kvar;
        li_tmp = li_tmp + 1;
        idx_tmp = iLoad_tmp.Next;
    end
    load_sequence = mod(0:nLoad-1, 24) + 1;
    DSSObj_tmp.delete;
    clear DSSObj_tmp DSSCircuit_tmp iLoad_tmp idx_tmp li_tmp;
    fprintf('[缓存] 负荷基准信息已重建: nLoad=%d\n', nLoad);
else
    fprintf('\n[计算] 缓存未命中，开始蒙特卡洛模拟...\n');

    % 初始化 OpenDSS COM 对象（仅初始化一次）
    DSSObj = actxserver('OpenDSSEngine.DSS');
    assert(DSSObj.Start(0), '[OpenDSS] 启动失败');
    DSSText    = DSSObj.Text;

    % 先编译一次获取系统信息
    DSSText.Command = ['Compile "', dssFile, '"'];
    DSSCircuit = DSSObj.ActiveCircuit;

    % 获取负荷列表
    iLoad = DSSCircuit.Loads;
    nLoad = iLoad.Count;
    loadNames = cell(nLoad, 1);
    loadKW_base  = zeros(nLoad, 1);
    loadKvar_base = zeros(nLoad, 1);
    idx = iLoad.First;
    li = 1;
    while idx > 0
        loadNames{li}    = iLoad.Name;
        loadKW_base(li)  = iLoad.kW;
        loadKvar_base(li) = iLoad.kvar;
        li = li + 1;
        idx = iLoad.Next;
    end
    fprintf('[OpenDSS] 检测到 %d 个负荷\n', nLoad);

    % 获取节点信息（用于确定nc）
    nodeNames = DSSCircuit.AllNodeNames;
    nc = length(nodeNames);
    busNames  = DSSCircuit.AllBusNames;

    % 获取线路信息
    iLine = DSSCircuit.Lines;
    nLine = iLine.Count;
    lineNames = cell(nLine, 1);
    idx = iLine.First;
    li = 1;
    while idx > 0
        lineNames{li} = iLine.Name;
        li = li + 1;
        idx = iLine.Next;
    end

    % 获取系统导纳矩阵（仅需一次，拓扑不变）
    sysY = DSSCircuit.SystemY;
    Y_vec = sysY(1:2:end) + 1i * sysY(2:2:end);
    nY = round(sqrt(length(Y_vec)));
    Y3ph = reshape(Y_vec, [nY, nY]);
    nodeOrderY = DSSCircuit.YNodeOrder;

    fprintf('[OpenDSS] 计算节点数 nc=%d, 线路数=%d\n', nc, nLine);
    fprintf('[OpenDSS] 导纳矩阵维度: %d x %d\n', nY, nY);

    % 负荷缩放序列：循环填充24种负荷
    load_sequence = mod(0:nLoad-1, 24) + 1;

    % 预分配存储
    V_mc_all     = zeros(N, nc);      % 复电压 (实部+j虚部)
    linePower_mc = zeros(N, nLine, 6); % 线路功率

    tic;
    h = waitbar(0, sprintf('三相MC潮流：0/%d (0.0%%)', N), 'Name', 'OpenDSS三相潮流');
    for t = 1:N
        % 重新编译电路（保证每次从基准状态开始）
        DSSText.Command = ['Compile "', dssFile, '"'];
        DSSCircuit = DSSObj.ActiveCircuit;

        % 设置各负荷的功率（按负荷曲线缩放）
        iLoad = DSSCircuit.Loads;
        idx = iLoad.First;
        li = 1;
        while idx > 0 && li <= nLoad
            scaleFactor = load_data{t, load_sequence(li)};
            iLoad.kW   = loadKW_base(li) * scaleFactor;
            iLoad.kvar = loadKvar_base(li) * scaleFactor;
            li = li + 1;
            idx = iLoad.Next;
        end

        % 求解潮流
        DSSSolution = DSSCircuit.Solution;
        DSSSolution.Solve;

        % 提取复电压
        allVolts = DSSCircuit.AllBusVolts;
        V_mc_all(t, :) = allVolts(1:2:end) + 1i * allVolts(2:2:end);

        % 提取线路功率
        iLineObj = DSSCircuit.Lines;
        idx2 = iLineObj.First;
        li2 = 1;
        while idx2 > 0
            DSSCircuit.SetActiveElement(['Line.' iLineObj.Name]);
            elem = DSSCircuit.ActiveCktElement;
            powers = elem.Powers;
            nTermPhases = length(powers) / 2;
            P_from = sum(powers(1:2:nTermPhases));
            Q_from = sum(powers(2:2:nTermPhases));
            P_to   = sum(powers(nTermPhases+1:2:end));
            Q_to   = sum(powers(nTermPhases+2:2:end));
            linePower_mc(t, li2, :) = [P_from, Q_from, P_to, Q_to, ...
                                        P_from+P_to, Q_from+Q_to];
            li2 = li2 + 1;
            idx2 = iLineObj.Next;
        end

        if mod(t, 100) == 0
            waitbar(t/N, h, sprintf('三相MC潮流：%d/%d (%.1f%%)', t, N, t/N*100));
        end
    end
    close(h);
    t1 = toc;

    % 保存缓存
    fprintf('[缓存] 保存结果到: %s\n', cacheFile);
    save(cacheFile, 'V_mc_all', 'nodeNames', 'nc', 'Y3ph', 'nodeOrderY', ...
         'loadNames', 'lineNames', 'linePower_mc', 'nLine', 'busNames', '-v7.3');
    fprintf('[缓存] 保存完成\n');
end

if t1 > 0
    fprintf('[MC] 蒙特卡洛模拟耗时: %.2f s\n', t1);
else
    fprintf('[MC] 已从缓存加载蒙特卡洛结果\n');
end
fprintf('[MC] 样本数 N=%d, 计算节点数 nc=%d\n', N, nc);

%% ========================================================================
%  段落2b: 全局绘图设置与节点编号
%  ========================================================================
% 全局绘图参数
fontSZ = 22;  fontSZ_title = 24;
fontName_cn = '宋体';  fontName_en = 'Times New Roman';
phaseChar_arr = {'A','B','C'};
colors_ph = {'b', 'r', [0 0.6 0]};  % A=蓝, B=红, C=绿

% 语义化节点标签: 从 nodeNames ('busname.phase') 生成 'BusName-Phase'
busRenameMap = containers.Map({'sourcebus','rg60'}, {'PCC','Reg'});
nodeLabels = cell(1, nc);
busLabels_of_node = cell(1, nc);   % 各节点对应的母线名（重命名后）
phaseLabels_of_node = zeros(1, nc); % 各节点对应的相序号 (1/2/3)
for i = 1:nc
    parts = split(nodeNames{i}, '.');
    busRaw = lower(parts{1});
    phNum  = str2double(parts{2});
    if busRenameMap.isKey(busRaw)
        busLabel = busRenameMap(busRaw);
    else
        busLabel = upper(busRaw);
    end
    busLabels_of_node{i}    = busLabel;
    phaseLabels_of_node(i)  = phNum;
    nodeLabels{i} = sprintf('%s-%s', busLabel, phaseChar_arr{phNum});
end

%% ========================================================================
%  段落3: 电压实同构（6n维 -> 主要内容式1-2）
%  ========================================================================
fprintf('\n[实同构] 构建 %d 维实电压向量...\n', 2*nc);

% V_mc_all: [N x nc] 复数矩阵
% 实同构: v_vec = [Re(V); Im(V)]，每行为一个样本的 2*nc 维实向量
v_vec = [real(V_mc_all), imag(V_mc_all)];  % [N x 2*nc]

% 各节点基准电压（用于标幺值显示）
V_base_node = mean(abs(V_mc_all), 1);  % [1 x nc] 各节点平均电压幅值
fprintf('[标幺] 各节点基准电压范围: %.1f ~ %.1f V\n', min(V_base_node), max(V_base_node));

% 电压统计量
V_mu    = mean(v_vec);        % 均值向量 [1 x 2*nc]
V_sigma = cov(v_vec);         % 协方差矩阵 [2*nc x 2*nc]
V_var   = diag(V_sigma);      % 方差向量

fprintf('[实同构] 均值向量维度: [1 x %d]\n', length(V_mu));
fprintf('[实同构] 协方差矩阵维度: [%d x %d]\n', size(V_sigma));

%% ========================================================================
%  段落4: GMM拟合（主要内容式10）
%  ========================================================================
fprintf('\n[GMM] 拟合 %d 分量的高斯混合模型...\n', gmm_num);

options = statset('Display', 'final', 'MaxIter', 500);
gmm_fit = gmdistribution.fit(v_vec, gmm_num, 'Regularize', 1e-10, ...
    'Options', options);
com_pi = gmm_fit.ComponentProportion;  % 各分量权重 [1 x K]

fprintf('[GMM] 拟合完成, BIC=%.2f\n', gmm_fit.BIC);

% Cholesky分解各分量协方差
B_GMM = cell(gmm_num, 1);
for k = 1:gmm_num
    B_GMM{k} = cholcov(gmm_fit.Sigma(:,:,k))';
end

% 从GMM采样用于对比
V_sample_GMM = [];
for k = 1:gmm_num
    nk = round(1e5 * com_pi(k));
    V_sample_GMM = [V_sample_GMM, ...
        mvnrnd(gmm_fit.mu(k,:)', gmm_fit.Sigma(:,:,k)+eps*eye(2*nc), nk)'];
end

fprintf('[GMM] GMM电压采样完成, 总样本数=%d\n', size(V_sample_GMM,2));

% --- Fig.2: 节点684电压实部GMM拟合（A/C两相，各一张独立图）---
fprintf('\n[图2] 绘制节点684电压实部GMM拟合...\n');
target_bus_fig2 = '684';
target_nodes_fig2 = find(cellfun(@(x) strcmp(x, upper(target_bus_fig2)), busLabels_of_node));
if ~isempty(target_nodes_fig2)
    for ii = 1:length(target_nodes_fig2)
        ni = target_nodes_fig2(ii);
        figure('Units','pixels','Position',[100+60*(ii-1),100+60*(ii-1),560,420]);
        Vb = V_base_node(ni);
        V_re_pu = v_vec(:, ni) / Vb;
        histogram(V_re_pu, 80, 'Normalization', 'pdf', 'FaceAlpha', 0.3, ...
            'FaceColor', [0.07 0.44 0.75], 'EdgeColor', 'none');
        hold on;
        x_plot_pu = linspace(min(V_re_pu)*0.999, max(V_re_pu)*1.001, 300);
        pdf_gmm_pu = zeros(size(x_plot_pu));
        for k = 1:gmm_num
            mu_k_pu = gmm_fit.mu(k, ni) / Vb;
            sig_k_pu = sqrt(gmm_fit.Sigma(ni, ni, k)) / Vb;
            pdf_gmm_pu = pdf_gmm_pu + com_pi(k) * normpdf(x_plot_pu, mu_k_pu, sig_k_pu);
        end
        plot(x_plot_pu, pdf_gmm_pu, 'r-', 'LineWidth', 4);
        hold off;
        set(gca,'fontsize',fontSZ,'fontname',fontName_en,'box','off');
        xlabel(sprintf('%s 电压实部(p.u.)', nodeLabels{ni}),'FontName',fontName_cn,'FontSize',fontSZ);
        ylabel('PDF','FontSize',fontSZ);
        legend({'MCS','GMM拟合'},'FontName',fontName_cn,'FontSize',fontSZ,'Location','best');
        saveas(gcf, fullfile(resultDir, sprintf('Fig02%c_VoltageRealGMM_%s.png', ...
            char('a'+ii-1), nodeLabels{ni})));
    end
    fprintf('[图2] 已保存\n');
else
    fprintf('[图2] 未找到母线 %s\n', target_bus_fig2);
end

%% ========================================================================
%  段落5: 三相二次型参数矩阵（主要内容式3-6）
%  ========================================================================
fprintf('\n[二次型] 构建三相功率二次型参数矩阵...\n');

% 构建节点-相映射（确保Y矩阵与电压向量的节点顺序一致）
% OpenDSS的 AllNodeNames 和 YNodeOrder 可能顺序不同，需要建立映射
% 将 YNodeOrder 转换为与 AllNodeNames 一致的排列
nodeOrderY_lower = lower(nodeOrderY);
nodeNames_lower  = lower(nodeNames);

% 建立 YNodeOrder -> AllNodeNames 的索引映射
permY2V = zeros(nc, 1);
for i = 1:nc
    idx_match = find(strcmp(nodeOrderY_lower{i}, nodeNames_lower));
    assert(~isempty(idx_match), '节点 %s 在 AllNodeNames 中未找到', nodeOrderY{i});
    permY2V(i) = idx_match;
end

% 将Y矩阵重排到与AllNodeNames一致的顺序
% permV2Y: V-node a 对应 Y-node permV2Y(a)
permV2Y = zeros(nc, 1);
permV2Y(permY2V) = 1:nc;
Y3ph_ordered = Y3ph(permV2Y, permV2Y);

% 构建节点注入功率二次型矩阵
[Node_P_matrix, Node_Q_matrix] = get_3ph_node_power_matrix(Y3ph_ordered, nc);

fprintf('[二次型] 节点功率矩阵: %d 个 [%d x %d]\n', nc, 2*nc, 2*nc);

% =========================================================================
%  构建线路连接信息（排除变压器连接）
% =========================================================================
% 从Y3ph_ordered中提取所有非零非对角的连接对
branchPairs_all = [];
for i = 1:nc
    for j = i+1:nc
        if abs(Y3ph_ordered(i,j)) > 1e-10
            branchPairs_all = [branchPairs_all; i, j];
        end
    end
end

% 排除变压器和电压调节器连接
% 判据1: 电压等级不同（如4.16kV vs 0.48kV），比值远大于1
% 判据2: 导纳幅值异常大（电压调节器虽然电压比≈1:1，但阻抗极低、导纳极大）
%         典型配电线路 |Y| < 50 S，调节器/短路连接 |Y| >> 50 S
Y_adm_threshold = 50;  % 导纳幅值阈值 (S)
isExcluded = false(size(branchPairs_all, 1), 1);
nXfm = 0; nReg = 0;
for l = 1:size(branchPairs_all, 1)
    ni = branchPairs_all(l,1); nj = branchPairs_all(l,2);
    Vratio = V_base_node(ni) / V_base_node(nj);
    Ymag = abs(Y3ph_ordered(ni, nj));
    if Vratio > 1.5 || Vratio < 1/1.5
        isExcluded(l) = true; nXfm = nXfm + 1;
    elseif Ymag > Y_adm_threshold
        isExcluded(l) = true; nReg = nReg + 1;
        fprintf('  [排除] 高导纳连接: %s → %s, |Y|=%.1f S (疑似调节器)\n', ...
            nodeNames{ni}, nodeNames{nj}, Ymag);
    end
end
fprintf('[线路识别] Y矩阵连接对: %d, 变压器: %d, 调节器: %d, 保留线路: %d\n', ...
    size(branchPairs_all,1), nXfm, nReg, sum(~isExcluded));

% 仅保留实际配电线路
branchPairs = branchPairs_all(~isExcluded, :);
nBranch = size(branchPairs, 1);

% 打印保留线路的导纳范围（确认过滤效果）
Y_mags_kept = zeros(nBranch, 1);
for l = 1:nBranch
    Y_mags_kept(l) = abs(Y3ph_ordered(branchPairs(l,1), branchPairs(l,2)));
end
fprintf('[线路识别] 保留线路 |Y| 范围: %.3f ~ %.3f S\n', min(Y_mags_kept), max(Y_mags_kept));

% 构建线路功率和损耗二次型矩阵（自导纳项，符号已修正 v2.0）
% 注: 此处仅使用线路自导纳 y_pp = -Y_bus(from_p, to_p)，
%     未包含相间互耦 y_pq (p≠q)。对于典型配电线路，互耦修正约10-30%。
[Line_P_matrix, Line_Q_matrix] = get_3ph_line_power_matrix(Y3ph_ordered, branchPairs, nc);
Line_loss_matrix = get_3ph_line_loss_matrix(Y3ph_ordered, branchPairs, nc);

fprintf('[二次型] 线路功率矩阵: %d 个 [%d x %d]\n', nBranch, 2*nc, 2*nc);

%% ========================================================================
%  段落6: GNCS参数计算（主要内容式14-21）
%  ========================================================================
fprintf('\n[GNCS] 计算线路功率的广义非中心卡方分布参数...\n');

% 线路GNCS参数
Line_P_lumda = cell(nBranch, gmm_num);
Line_P_delta = cell(nBranch, gmm_num);
Line_P_const = cell(nBranch, gmm_num);
Line_Q_lumda = cell(nBranch, gmm_num);
Line_Q_delta = cell(nBranch, gmm_num);
Line_Q_const = cell(nBranch, gmm_num);

tic;
h = waitbar(0, sprintf('线路GNCS参数：0/%d', nBranch), 'Name', '线路参数计算');
for l = 1:nBranch
    for k = 1:gmm_num
        mu_k = gmm_fit.mu(k,:)';

        % 有功
        BMB = B_GMM{k}' * Line_P_matrix{l} * B_GMM{k};
        BMB_sym = 0.5*(BMB + BMB');
        [D_mat, Lambda_mat] = eig(BMB_sym);
        lv = diag(Lambda_mat);
        a_P = 2 * B_GMM{k}' * Line_P_matrix{l} * mu_k;
        bv = D_mat' * a_P;
        alpha_P = mu_k' * Line_P_matrix{l} * mu_k;
        tol = 1e-12 * max(abs(lv));
        inz = abs(lv) > tol;
        cv = zeros(size(lv));
        cv(inz) = bv(inz) ./ (2*lv(inz));
        dv = cv.^2;
        Line_P_lumda{l,k} = lv;
        Line_P_delta{l,k} = dv;
        Line_P_const{l,k} = alpha_P - sum(lv(inz).*dv(inz));

        % 无功
        BMB_Q = B_GMM{k}' * Line_Q_matrix{l} * B_GMM{k};
        BMB_Q_sym = 0.5*(BMB_Q + BMB_Q');
        [D_Q, L_Q] = eig(BMB_Q_sym);
        lq = diag(L_Q);
        a_Q = 2 * B_GMM{k}' * Line_Q_matrix{l} * mu_k;
        bq = D_Q' * a_Q;
        alpha_Q = mu_k' * Line_Q_matrix{l} * mu_k;
        tol_q = 1e-12 * max(abs(lq));
        inz_q = abs(lq) > tol_q;
        cq = zeros(size(lq));
        cq(inz_q) = bq(inz_q) ./ (2*lq(inz_q));
        dq = cq.^2;
        Line_Q_lumda{l,k} = lq;
        Line_Q_delta{l,k} = dq;
        Line_Q_const{l,k} = alpha_Q - sum(lq(inz_q).*dq(inz_q));
    end
    waitbar(l/nBranch, h, sprintf('线路GNCS参数：%d/%d', l, nBranch));
end
close(h);
t3 = toc;
fprintf('[GNCS] 线路参数计算耗时: %.2f s\n', t3);

%% ========================================================================
%  段落7: 线路功率统计矩计算与误差分析（主要内容式24-29）
%  ========================================================================
fprintf('\n[统计矩] 计算线路功率均值/方差并与GMM采样对比...\n');

% 筛选同相线路连接（排除跨相互耦项，仅保留 phase_from == phase_to）
samePhaseIdx = [];
for l = 1:nBranch
    parts_i = split(nodeNames{branchPairs(l,1)}, '.');
    parts_j = split(nodeNames{branchPairs(l,2)}, '.');
    if strcmp(parts_i{2}, parts_j{2})
        samePhaseIdx = [samePhaseIdx; l];
    end
end
nSP = length(samePhaseIdx);
fprintf('[统计矩] 同相线路连接数: %d / %d\n', nSP, nBranch);

% MC采样直接计算线路功率（作为基准，使用修正后的二次型矩阵）
V_mc_T = v_vec';  % [2*nc x N]，MC电压实同构转置
Line_Power_P_MC = cell(nBranch, 1);
Line_Power_Q_MC = cell(nBranch, 1);
for l = 1:nBranch
    Line_Power_P_MC{l} = dot(V_mc_T, Line_P_matrix{l} * V_mc_T)';
    Line_Power_Q_MC{l} = dot(V_mc_T, Line_Q_matrix{l} * V_mc_T)';
end

% 校验1: 打印线路功率范围确认物理合理性
allP_mean = cellfun(@mean, Line_Power_P_MC);
allQ_mean = cellfun(@mean, Line_Power_Q_MC);
fprintf('[校验] 二次型法线路有功范围: %.2f ~ %.2f kW\n', min(allP_mean)/1e3, max(allP_mean)/1e3);
fprintf('[校验] 二次型法线路无功范围: %.2f ~ %.2f kVar\n', min(allQ_mean)/1e3, max(allQ_mean)/1e3);

% 校验2: 用复电压直接计算线路功率，与二次型结果对比
% P_ij = Re(V_i * conj(y_ij * (V_i - V_j))), y_ij = -Y_bus(i,j)
fprintf('[校验] 复电压直接计算线路功率进行交叉验证...\n');
nCheckLines = min(5, nBranch);  % 检查前几条线路
for l = 1:nCheckLines
    ni = branchPairs(l,1); nj = branchPairs(l,2);
    y_ij = -Y3ph_ordered(ni, nj);  % 线路原始导纳
    Vi = V_mc_all(:, ni);  % [N x 1] 复电压
    Vj = V_mc_all(:, nj);
    S_direct = Vi .* conj(y_ij .* (Vi - Vj));  % [N x 1]
    P_direct_mean = mean(real(S_direct));
    P_qform_mean  = mean(Line_Power_P_MC{l});
    fprintf('  线路%d (%s→%s): 直接法=%.2f kW, 二次型=%.2f kW, 差异=%.4f%%\n', ...
        l, nodeNames{ni}, nodeNames{nj}, P_direct_mean/1e3, P_qform_mean/1e3, ...
        abs(P_direct_mean - P_qform_mean)/max(abs(P_direct_mean), 1)*100);
end

% 解析计算各分量的均值和方差（主要内容式24-25）
mean_Line_P_com = zeros(nBranch, gmm_num);
mean_Line_Q_com = zeros(nBranch, gmm_num);
var_Line_P_com  = zeros(nBranch, gmm_num);
var_Line_Q_com  = zeros(nBranch, gmm_num);

for l = 1:nBranch
    for k = 1:gmm_num
        % 均值: E[P] = tr(M*Sigma) + mu'*M*mu（式24）
        mean_Line_P_com(l,k) = trace(Line_P_matrix{l} * gmm_fit.Sigma(:,:,k)) ...
            + gmm_fit.mu(k,:) * Line_P_matrix{l} * gmm_fit.mu(k,:)';
        mean_Line_Q_com(l,k) = trace(Line_Q_matrix{l} * gmm_fit.Sigma(:,:,k)) ...
            + gmm_fit.mu(k,:) * Line_Q_matrix{l} * gmm_fit.mu(k,:)';

        % 方差: Var = 2*tr((M*Sigma)^2) + 4*mu'*M*Sigma*M*mu（式25）
        var_Line_P_com(l,k) = 2*trace((Line_P_matrix{l}*gmm_fit.Sigma(:,:,k))^2) ...
            + 4*gmm_fit.mu(k,:)*Line_P_matrix{l}*gmm_fit.Sigma(:,:,k)*Line_P_matrix{l}*gmm_fit.mu(k,:)';
        var_Line_Q_com(l,k) = 2*trace((Line_Q_matrix{l}*gmm_fit.Sigma(:,:,k))^2) ...
            + 4*gmm_fit.mu(k,:)*Line_Q_matrix{l}*gmm_fit.Sigma(:,:,k)*Line_Q_matrix{l}*gmm_fit.mu(k,:)';
    end
end

% GMM混合后的均值和方差（式26-27）
mean_Line_P_this = mean_Line_P_com * com_pi';
mean_Line_Q_this = mean_Line_Q_com * com_pi';
var_Line_P_this  = var_Line_P_com * com_pi' + ...
    sum(com_pi .* (mean_Line_P_com - mean_Line_P_this).^2, 2);
var_Line_Q_this  = var_Line_Q_com * com_pi' + ...
    sum(com_pi .* (mean_Line_Q_com - mean_Line_Q_this).^2, 2);

% MC采样基准
mean_Line_P_mc = cellfun(@mean, Line_Power_P_MC);
mean_Line_Q_mc = cellfun(@mean, Line_Power_Q_MC);
var_Line_P_mc  = cellfun(@var,  Line_Power_P_MC);
var_Line_Q_mc  = cellfun(@var,  Line_Power_Q_MC);

% 相对误差（仅对同相连接计算）
% 使用 max(|mean|, threshold) 作为分母，避免小分母放大效应
% threshold = 各量绝对值中位数的1%，确保分母有意义
mean_P_sp = mean_Line_P_mc(samePhaseIdx);
mean_Q_sp = mean_Line_Q_mc(samePhaseIdx);
var_P_sp  = var_Line_P_mc(samePhaseIdx);
var_Q_sp  = var_Line_Q_mc(samePhaseIdx);

thr_P = 0.01 * median(abs(mean_P_sp));
thr_Q = 0.01 * median(abs(mean_Q_sp));
thr_varP = 0.01 * median(abs(var_P_sp));
thr_varQ = 0.01 * median(abs(var_Q_sp));

error_Line_P_mean = abs(mean_P_sp - mean_Line_P_this(samePhaseIdx)) ...
    ./ max(abs(mean_P_sp), thr_P) * 100;
error_Line_Q_mean = abs(mean_Q_sp - mean_Line_Q_this(samePhaseIdx)) ...
    ./ max(abs(mean_Q_sp), thr_Q) * 100;
error_Line_P_var  = abs(var_P_sp  - var_Line_P_this(samePhaseIdx)) ...
    ./ max(abs(var_P_sp), thr_varP) * 100;
error_Line_Q_var  = abs(var_Q_sp  - var_Line_Q_this(samePhaseIdx)) ...
    ./ max(abs(var_Q_sp), thr_varQ) * 100;

% 处理可能残留的 NaN/Inf
error_Line_P_mean(~isfinite(error_Line_P_mean)) = 0;
error_Line_Q_mean(~isfinite(error_Line_Q_mean)) = 0;
error_Line_P_var(~isfinite(error_Line_P_var))   = 0;
error_Line_Q_var(~isfinite(error_Line_Q_var))   = 0;

fprintf('[统计矩] 线路有功均值最大误差: %.4f%%\n', max(error_Line_P_mean));
fprintf('[统计矩] 线路无功均值最大误差: %.4f%%\n', max(error_Line_Q_mean));
fprintf('[统计矩] 线路有功方差最大误差: %.4f%%\n', max(error_Line_P_var));
fprintf('[统计矩] 线路无功方差最大误差: %.4f%%\n', max(error_Line_Q_var));

% 输出绝对误差供参考（识别小分母线路）
absErr_P_mean = abs(mean_P_sp - mean_Line_P_this(samePhaseIdx));
absErr_Q_mean = abs(mean_Q_sp - mean_Line_Q_this(samePhaseIdx));
[~, worstP] = max(error_Line_P_mean);
[~, worstQ] = max(error_Line_Q_mean);
fprintf('[统计矩] 有功均值最大相对误差线路: 均值=%.4e, 绝对误差=%.4e\n', ...
    mean_P_sp(worstP), absErr_P_mean(worstP));
fprintf('[统计矩] 无功均值最大相对误差线路: 均值=%.4e, 绝对误差=%.4e\n', ...
    mean_Q_sp(worstQ), absErr_Q_mean(worstQ));

% 构建线路标签（使用语义化节点标签）
branchLabels = cell(nSP, 1);
for idx = 1:nSP
    l = samePhaseIdx(idx);
    ni = branchPairs(l,1);  nj = branchPairs(l,2);
    branchLabels{idx} = sprintf('%s→%s', nodeLabels{ni}, nodeLabels{nj});
end

% (旧 Fig.6a/6b 线路功率误差柱状图已删除，误差汇总见段落12表1)

% --- 三相节点识别（用于Fig.10和段落10）---
busPhaseMap = containers.Map();
for i = 1:nc
    parts = split(nodeNames{i}, '.');
    busname = parts{1};
    phase = str2double(parts{2});
    if ~busPhaseMap.isKey(busname)
        busPhaseMap(busname) = [0, 0, 0];
    end
    phaseIdx_tmp = busPhaseMap(busname);
    phaseIdx_tmp(phase) = i;
    busPhaseMap(busname) = phaseIdx_tmp;
end
allBuses = keys(busPhaseMap);
threePhaseBuses = {}; threePhaseBusIdx = [];
for b = 1:length(allBuses)
    phaseIdx_tmp = busPhaseMap(allBuses{b});
    if all(phaseIdx_tmp > 0)
        threePhaseBuses{end+1} = allBuses{b};
        threePhaseBusIdx = [threePhaseBusIdx; phaseIdx_tmp];
    end
end
n3ph = length(threePhaseBuses);
fprintf('[三相识别] 三相完整节点数: %d\n', n3ph);

%% ========================================================================
%  段落8: 线路功率概率分布（主要内容式14/21/35, Davies法）
%  ========================================================================
fprintf('\n[PDF] 计算线路功率概率分布...\n');

% 按物理线路（母线对）分组同相连接，识别三相/两相/单相线路
physLineMap = containers.Map();
for idx = 1:nSP
    l = samePhaseIdx(idx);
    ni = branchPairs(l,1); nj = branchPairs(l,2);
    pi_parts = split(nodeNames{ni}, '.');
    pj_parts = split(nodeNames{nj}, '.');
    % 用母线名作为key（同一物理线路不同相的母线名相同）
    busFrom = lower(pi_parts{1}); busTo = lower(pj_parts{1});
    % 按字典序排列确保唯一key
    if string(busFrom) > string(busTo), lineKey = [busTo '-' busFrom];
    else, lineKey = [busFrom '-' busTo]; end
    ph = str2double(pi_parts{2});
    if ~physLineMap.isKey(lineKey)
        physLineMap(lineKey) = struct('brIdx', [0,0,0], ...
            'from_nodes', [0,0,0], 'to_nodes', [0,0,0], ...
            'from_bus', pi_parts{1}, 'to_bus', pj_parts{1});
    end
    info = physLineMap(lineKey);
    info.brIdx(ph) = l;
    info.from_nodes(ph) = ni; info.to_nodes(ph) = nj;
    physLineMap(lineKey) = info;
end

lineKeys_all = keys(physLineMap);
threePhaseLineKeys = {}; twoPhaseLineKeys = {}; singlePhaseLineKeys = {};
for idx = 1:length(lineKeys_all)
    info = physLineMap(lineKeys_all{idx});
    nPh = sum(info.brIdx > 0);
    if nPh == 3,     threePhaseLineKeys{end+1}  = lineKeys_all{idx};
    elseif nPh == 2, twoPhaseLineKeys{end+1}    = lineKeys_all{idx};
    elseif nPh == 1, singlePhaseLineKeys{end+1} = lineKeys_all{idx};
    end
end
fprintf('[PDF] 三相: %d, 两相: %d, 单相: %d\n', ...
    length(threePhaseLineKeys), length(twoPhaseLineKeys), length(singlePhaseLineKeys));

n_pdf_pts = 200;
n_samp_fig = 1e5;

% (旧 Fig.7/8/9/10 已删除，线路功率PDF对比见段落12新图4)


%% ========================================================================
%  段落9: 三相电压幅值分析（主要内容式38-46b + 均值方差精确/近似）
%  ========================================================================
fprintf('\n[电压幅值] 计算三相电压幅值概率分布与均值方差...\n');

% 9.1 构造电压幅值平方的二次型矩阵 |V_i|^2 = V_e^2 + V_f^2
M_Vmag2 = cell(nc, 1);
for i = 1:nc
    M_Vmag2{i} = zeros(2*nc, 2*nc);
    M_Vmag2{i}(i, i) = 1;         % 实部^2
    M_Vmag2{i}(nc+i, nc+i) = 1;   % 虚部^2
end

% 9.2 计算 |V|^2 的GNCS参数
Vmag2_lumda = cell(nc, gmm_num);
Vmag2_delta = cell(nc, gmm_num);
Vmag2_const = cell(nc, gmm_num);

for i = 1:nc
    for k = 1:gmm_num
        BMB = B_GMM{k}' * M_Vmag2{i} * B_GMM{k};
        BMB_sym = 0.5*(BMB + BMB');
        [Q_v, D_v] = eig(BMB_sym);
        lv = diag(D_v);
        mu_k = gmm_fit.mu(k,:)';
        a_v = 2 * B_GMM{k}' * M_Vmag2{i} * mu_k;
        bv = Q_v' * a_v;
        alpha_v = mu_k' * M_Vmag2{i} * mu_k;

        tol = 1e-6 * max(abs(lv));
        inz = abs(lv) > tol;
        cv = zeros(size(lv));
        cv(inz) = bv(inz) ./ (2*lv(inz));
        dv = cv.^2;

        % 仅保留有效特征值（M_Vmag2 rank=2，过滤近零项避免数值问题）
        Vmag2_lumda{i,k} = lv(inz);
        Vmag2_delta{i,k} = dv(inz);
        Vmag2_const{i,k} = alpha_v - sum(lv(inz).*dv(inz));
    end
end

% 9.3 MC基准: 电压幅值
Vmag_mc = zeros(N, nc);
for i = 1:nc
    Vmag_mc(:,i) = sqrt(v_vec(:,i).^2 + v_vec(:,nc+i).^2);
end
mean_Vmag_mc = mean(Vmag_mc, 1)';
var_Vmag_mc  = var(Vmag_mc, 0, 1)';

% 9.4 精确计算: E[|V|] 和 Var(|V|) 通过GNCS配方采样
%     原理：从 |V|^2 的 GNCS 分布直接采样 -> 取 sqrt 得 |V| 样本 -> 计算矩
%     避免 Davies 法在大数值(|V|^2~10^6)下的振荡积分失准问题
n_gncs_samples = 1e5;  % GNCS采样数（足够精确且快速）
mean_Vmag_exact = zeros(nc, 1);
var_Vmag_exact  = zeros(nc, 1);

tic_vmag_exact = tic;
for i = 1:nc
    Vmag_gncs_samples = [];
    for k = 1:gmm_num
        nk = round(n_gncs_samples * com_pi(k));
        if nk == 0, continue; end
        n_eig = length(Vmag2_lumda{i,k});
        Z = randn(n_eig, nk);
        c_vals = sqrt(Vmag2_delta{i,k});  % c = sqrt(delta)
        % |V|^2 样本 = const + sum(lambda_j * (z_j + c_j)^2)
        Vsq_samples = Vmag2_const{i,k} + sum(Vmag2_lumda{i,k} .* (Z + c_vals).^2, 1);
        % 取 sqrt 得 |V| 样本（确保非负）
        Vmag_gncs_samples = [Vmag_gncs_samples, sqrt(max(Vsq_samples, 0))];
    end
    mean_Vmag_exact(i) = mean(Vmag_gncs_samples);
    var_Vmag_exact(i)  = var(Vmag_gncs_samples);
end
t_vmag_exact = toc(tic_vmag_exact);
fprintf('[电压幅值] GNCS配方采样完成 (每节点%d样本), 耗时: %.2f s\n', n_gncs_samples, t_vmag_exact);

% 9.5 近似计算: Delta法
% E[|V|] ≈ sqrt(E[|V|^2]),  Var(|V|) ≈ Var(|V|^2) / (4*E[|V|^2])
mean_Vmag2_gncs = zeros(nc, 1);
var_Vmag2_gncs  = zeros(nc, 1);

tic_vmag_delta = tic;
for i = 1:nc
    mean_W_k = zeros(gmm_num, 1);
    var_W_k  = zeros(gmm_num, 1);
    for k = 1:gmm_num
        mean_W_k(k) = Vmag2_const{i,k} + sum(Vmag2_lumda{i,k} .* (1 + Vmag2_delta{i,k}));
        var_W_k(k)  = 2 * sum(Vmag2_lumda{i,k}.^2 .* (1 + 2*Vmag2_delta{i,k}));
    end
    mean_Vmag2_gncs(i) = mean_W_k' * com_pi';
    var_Vmag2_gncs(i)  = var_W_k' * com_pi' + mean_W_k.^2' * com_pi' - mean_Vmag2_gncs(i)^2;
end

mean_Vmag_approx = sqrt(mean_Vmag2_gncs);
var_Vmag_approx  = var_Vmag2_gncs ./ (4 * mean_Vmag2_gncs);
t_vmag_delta = toc(tic_vmag_delta);
fprintf('[电压幅值] Delta法完成, 耗时: %.4f s\n', t_vmag_delta);

% 9.6 误差对比
error_Vmag_mean_exact  = abs((mean_Vmag_mc - mean_Vmag_exact) ./ mean_Vmag_mc) * 100;
error_Vmag_mean_approx = abs((mean_Vmag_mc - mean_Vmag_approx) ./ mean_Vmag_mc) * 100;
error_Vmag_var_exact   = abs((var_Vmag_mc - var_Vmag_exact) ./ var_Vmag_mc) * 100;
error_Vmag_var_approx  = abs((var_Vmag_mc - var_Vmag_approx) ./ var_Vmag_mc) * 100;

fprintf('[电压幅值] 均值精确最大误差: %.4f%%, 近似最大误差: %.4f%%\n', ...
    max(error_Vmag_mean_exact), max(error_Vmag_mean_approx));
fprintf('[电压幅值] 方差精确最大误差: %.4f%%, 近似最大误差: %.4f%%\n', ...
    max(error_Vmag_var_exact), max(error_Vmag_var_approx));

% (Fig.4a/4b/4c/4d 电压幅值误差图在段落12与SOTA对比后统一绘制)
% (Fig.5 节点26-28电压幅值PDF已移除，电压幅值PDF对比见段落12的Fig.4e/4f/4g)

%% ========================================================================
%  段落10: VUF分析（主要内容式47-56 + 均值方差精确/近似）
%  ========================================================================
fprintf('\n[VUF] 计算电压不平衡度的概率分布与均值方差...\n');
% 10.1 三相节点已在段落8前识别 (threePhaseBuses, threePhaseBusIdx, n3ph)
fprintf('[VUF] 使用已识别的 %d 个三相节点\n', n3ph);

% 10.2 对称分量变换矩阵
a = exp(1i * 2*pi/3);  % 120度相移因子
T_seq = (1/3) * [1, 1, 1; 1, a, a^2; 1, a^2, a];  % 正/负/零序变换

% 10.3 MC基准: 计算每个三相节点的VUF
VUF_mc = zeros(N, n3ph);
Vmag_pos_mc = zeros(N, n3ph);
Vmag_neg_mc = zeros(N, n3ph);
Vmag_zero_mc = zeros(N, n3ph);

for b = 1:n3ph
    idxA = threePhaseBusIdx(b, 1);
    idxB = threePhaseBusIdx(b, 2);
    idxC = threePhaseBusIdx(b, 3);

    V_abc = [V_mc_all(:, idxA), V_mc_all(:, idxB), V_mc_all(:, idxC)];  % [N x 3]
    V_seq = (T_seq * V_abc.')';  % [N x 3]: 列1=零序V0, 列2=正序V+, 列3=负序V-

    Vmag_pos_mc(:, b) = abs(V_seq(:, 2));  % 正序幅值 V+
    Vmag_neg_mc(:, b) = abs(V_seq(:, 3));  % 负序幅值 V-
    Vmag_zero_mc(:, b) = abs(V_seq(:, 1)); % 零序幅值 V0

    VUF_mc(:, b) = Vmag_neg_mc(:, b) ./ Vmag_pos_mc(:, b) * 100;  % 百分比
end

mean_VUF_mc = mean(VUF_mc, 1)';
var_VUF_mc  = var(VUF_mc, 0, 1)';

% MC基准: 序电压幅值统计量
mean_Vpos_mc = mean(Vmag_pos_mc, 1)'; var_Vpos_mc = var(Vmag_pos_mc, 0, 1)';
mean_Vneg_mc = mean(Vmag_neg_mc, 1)'; var_Vneg_mc = var(Vmag_neg_mc, 0, 1)';
mean_V0_mc   = mean(Vmag_zero_mc, 1)'; var_V0_mc = var(Vmag_zero_mc, 0, 1)';

% 10.4 构造序电压幅值平方的二次型矩阵
% |V+|^2, |V-|^2, |V0|^2 都可以表示为实电压的二次型
M_Vpos2 = cell(n3ph, 1);
M_Vneg2 = cell(n3ph, 1);
M_V0sq  = cell(n3ph, 1);  % 零序

for b = 1:n3ph
    idxA = threePhaseBusIdx(b, 1);
    idxB = threePhaseBusIdx(b, 2);
    idxC = threePhaseBusIdx(b, 3);

    c_pos = [1, a, a^2] / 3;
    c_neg = [1, a^2, a] / 3;
    c_zero = [1, 1, 1] / 3;

    idx_e = [idxA, idxB, idxC];
    idx_f = [nc+idxA, nc+idxB, nc+idxC];

    cr_pos = real(c_pos); ci_pos = imag(c_pos);
    cr_neg = real(c_neg); ci_neg = imag(c_neg);
    cr_zero = real(c_zero); ci_zero = imag(c_zero);

    M_pos = zeros(2*nc, 2*nc);
    M_neg = zeros(2*nc, 2*nc);
    M_zero = zeros(2*nc, 2*nc);

    for p = 1:3
        for q = 1:3
            M_pos(idx_e(p), idx_e(q)) = cr_pos(p)*cr_pos(q) + ci_pos(p)*ci_pos(q);
            M_pos(idx_e(p), idx_f(q)) = ci_pos(p)*cr_pos(q) - cr_pos(p)*ci_pos(q);
            M_pos(idx_f(p), idx_e(q)) = cr_pos(p)*ci_pos(q) - ci_pos(p)*cr_pos(q);
            M_pos(idx_f(p), idx_f(q)) = cr_pos(p)*cr_pos(q) + ci_pos(p)*ci_pos(q);

            M_neg(idx_e(p), idx_e(q)) = cr_neg(p)*cr_neg(q) + ci_neg(p)*ci_neg(q);
            M_neg(idx_e(p), idx_f(q)) = ci_neg(p)*cr_neg(q) - cr_neg(p)*ci_neg(q);
            M_neg(idx_f(p), idx_e(q)) = cr_neg(p)*ci_neg(q) - ci_neg(p)*cr_neg(q);
            M_neg(idx_f(p), idx_f(q)) = cr_neg(p)*cr_neg(q) + ci_neg(p)*ci_neg(q);

            M_zero(idx_e(p), idx_e(q)) = cr_zero(p)*cr_zero(q) + ci_zero(p)*ci_zero(q);
            M_zero(idx_e(p), idx_f(q)) = ci_zero(p)*cr_zero(q) - cr_zero(p)*ci_zero(q);
            M_zero(idx_f(p), idx_e(q)) = cr_zero(p)*ci_zero(q) - ci_zero(p)*cr_zero(q);
            M_zero(idx_f(p), idx_f(q)) = cr_zero(p)*cr_zero(q) + ci_zero(p)*ci_zero(q);
        end
    end

    M_Vpos2{b} = 0.5 * (M_pos + M_pos');
    M_Vneg2{b} = 0.5 * (M_neg + M_neg');
    M_V0sq{b}  = 0.5 * (M_zero + M_zero');
end

% 10.5 计算序电压幅值的GNCS参数 (正序、负序、零序)
Vpos2_lumda = cell(n3ph, gmm_num);  Vpos2_delta = cell(n3ph, gmm_num);  Vpos2_const = cell(n3ph, gmm_num);
Vneg2_lumda = cell(n3ph, gmm_num);  Vneg2_delta = cell(n3ph, gmm_num);  Vneg2_const = cell(n3ph, gmm_num);
V0sq_lumda  = cell(n3ph, gmm_num);  V0sq_delta  = cell(n3ph, gmm_num);  V0sq_const  = cell(n3ph, gmm_num);

for b = 1:n3ph
    for k = 1:gmm_num
        mu_k = gmm_fit.mu(k,:)';

        % 正序（过滤近零特征值，M_Vpos2 rank<=6）
        BMB = B_GMM{k}' * M_Vpos2{b} * B_GMM{k};
        BMB_sym = 0.5*(BMB+BMB');
        [Q,D] = eig(BMB_sym); lv = diag(D);
        a_v = 2*B_GMM{k}'*M_Vpos2{b}*mu_k;
        bv = Q'*a_v; alpha_v = mu_k'*M_Vpos2{b}*mu_k;
        tol = 1e-6*max(abs(lv)); inz = abs(lv)>tol;
        cv = zeros(size(lv)); cv(inz) = bv(inz)./(2*lv(inz));
        dv = cv.^2;
        Vpos2_lumda{b,k} = lv(inz); Vpos2_delta{b,k} = dv(inz);
        Vpos2_const{b,k} = alpha_v - sum(lv(inz).*dv(inz));

        % 负序（过滤近零特征值，M_Vneg2 rank<=6）
        BMB = B_GMM{k}' * M_Vneg2{b} * B_GMM{k};
        BMB_sym = 0.5*(BMB+BMB');
        [Q,D] = eig(BMB_sym); lv = diag(D);
        a_v = 2*B_GMM{k}'*M_Vneg2{b}*mu_k;
        bv = Q'*a_v; alpha_v = mu_k'*M_Vneg2{b}*mu_k;
        tol = 1e-6*max(abs(lv)); inz = abs(lv)>tol;
        cv = zeros(size(lv)); cv(inz) = bv(inz)./(2*lv(inz));
        dv = cv.^2;
        Vneg2_lumda{b,k} = lv(inz); Vneg2_delta{b,k} = dv(inz);
        Vneg2_const{b,k} = alpha_v - sum(lv(inz).*dv(inz));

        % 零序（过滤近零特征值，M_V0sq rank<=6）
        BMB = B_GMM{k}' * M_V0sq{b} * B_GMM{k};
        BMB_sym = 0.5*(BMB+BMB');
        [Q,D] = eig(BMB_sym); lv = diag(D);
        a_v = 2*B_GMM{k}'*M_V0sq{b}*mu_k;
        bv = Q'*a_v; alpha_v = mu_k'*M_V0sq{b}*mu_k;
        tol = 1e-6*max(abs(lv)); inz = abs(lv)>tol;
        cv = zeros(size(lv)); cv(inz) = bv(inz)./(2*lv(inz));
        dv = cv.^2;
        V0sq_lumda{b,k} = lv(inz); V0sq_delta{b,k} = dv(inz);
        V0sq_const{b,k} = alpha_v - sum(lv(inz).*dv(inz));
    end
end

% 10.6 VUF解析PDF参数计算（主要内容.md 式47-51，附录A.2）
%     对每个GMM分量k和每个三相母线b，计算:
%       λ+: 标准化正序非中心参数 = |E[V+]|²/σ²+
%       λ-: 标准化负序非中心参数 = |E[V-]|²/σ²-
%       ρ:  Cor(|V+|², |V-|²)
%       σ_ratio: σ-/σ+ (尺度比)
fprintf('[VUF] 计算解析PDF参数（主要内容.md 式51）...\n');

VUF_lam_pos   = zeros(n3ph, gmm_num);
VUF_lam_neg   = zeros(n3ph, gmm_num);
VUF_rho       = zeros(n3ph, gmm_num);
VUF_sig_ratio = zeros(n3ph, gmm_num);

for b = 1:n3ph
    idxA = threePhaseBusIdx(b, 1);
    idxB = threePhaseBusIdx(b, 2);
    idxC = threePhaseBusIdx(b, 3);

    % 对称分量系数
    c_pos_coeff = [1, a, a^2] / 3;   % 正序
    c_neg_coeff = [1, a^2, a] / 3;   % 负序

    % 构造线性映射 L_pos, L_neg (2 x 2nc)
    %   [Re(V+), Im(V+)]' = L_pos * V_real
    idx_e = [idxA, idxB, idxC];
    idx_f = [nc+idxA, nc+idxB, nc+idxC];

    L_pos = zeros(2, 2*nc);
    L_pos(1, idx_e) = real(c_pos_coeff);
    L_pos(1, idx_f) = -imag(c_pos_coeff);
    L_pos(2, idx_e) = imag(c_pos_coeff);
    L_pos(2, idx_f) = real(c_pos_coeff);

    L_neg = zeros(2, 2*nc);
    L_neg(1, idx_e) = real(c_neg_coeff);
    L_neg(1, idx_f) = -imag(c_neg_coeff);
    L_neg(2, idx_e) = imag(c_neg_coeff);
    L_neg(2, idx_f) = real(c_neg_coeff);

    for k = 1:gmm_num
        mu_k = gmm_fit.mu(k,:)';
        Sigma_k = gmm_fit.Sigma(:,:,k);

        % 正序: 均值和协方差
        mu_pos = L_pos * mu_k;               % [Re(E[V+]), Im(E[V+])]
        Sigma_pos = L_pos * Sigma_k * L_pos'; % 2x2 协方差
        sig2_pos = trace(Sigma_pos) / 2;      % 平均特征值（近似圆对称）
        VUF_lam_pos(b,k) = (mu_pos' * mu_pos) / sig2_pos;

        % 负序: 均值和协方差
        mu_neg = L_neg * mu_k;
        Sigma_neg = L_neg * Sigma_k * L_neg';
        sig2_neg = trace(Sigma_neg) / 2;
        VUF_lam_neg(b,k) = (mu_neg' * mu_neg) / sig2_neg;

        % 尺度比
        VUF_sig_ratio(b,k) = sqrt(sig2_neg / sig2_pos);

        % 相关系数 ρ = Cor(|V+|², |V-|²)
        % 利用式(29): Cov(Q1,Q2) = 2·tr(M1·Σ·M2·Σ) + 4·μ'·M1·Σ·M2·μ
        MSp = M_Vpos2{b} * Sigma_k;
        MSn = M_Vneg2{b} * Sigma_k;
        cov_pn = 2*trace(MSp * MSn) + 4*mu_k'*M_Vpos2{b}*Sigma_k*M_Vneg2{b}*mu_k;
        var_p  = 2*trace(MSp^2)     + 4*mu_k'*M_Vpos2{b}*Sigma_k*M_Vpos2{b}*mu_k;
        var_n  = 2*trace(MSn^2)     + 4*mu_k'*M_Vneg2{b}*Sigma_k*M_Vneg2{b}*mu_k;

        rho_k = cov_pn / sqrt(var_p * var_n);
        VUF_rho(b,k) = max(min(rho_k, 0.999), -0.999);
    end
end
fprintf('[VUF] PDF参数计算完成: λ+, λ-, ρ, σ_ratio\n');

% 10.7 精确计算VUF均值方差（GMM电压采样，保留V+/V-相关性）
%     直接从GMM电压样本计算序电压，然后得到VUF
%     注意：不能独立采样|V+|²和|V-|²，因为来自同一个电压向量，存在相关性
%     使用5e5样本以减小方差估计的有限样本误差
n_vuf_exact_samp = 5e5;
fprintf('[VUF] GMM电压采样计算精确矩（N=%d，保留相关性）...\n', n_vuf_exact_samp);
mean_VUF_exact = zeros(n3ph, 1);
var_VUF_exact  = zeros(n3ph, 1);

% 生成大量GMM电压样本（专用于VUF精确矩计算）
V_sample_VUF = [];
for k = 1:gmm_num
    nk = round(n_vuf_exact_samp * com_pi(k));
    V_sample_VUF = [V_sample_VUF, ...
        mvnrnd(gmm_fit.mu(k,:)', gmm_fit.Sigma(:,:,k)+eps*eye(2*nc), nk)'];
end

for b = 1:n3ph
    idxA = threePhaseBusIdx(b, 1);
    idxB = threePhaseBusIdx(b, 2);
    idxC = threePhaseBusIdx(b, 3);

    % 从GMM电压样本提取三相复电压
    Va_gmm = V_sample_VUF(idxA, :) + 1i * V_sample_VUF(nc+idxA, :);
    Vb_gmm = V_sample_VUF(idxB, :) + 1i * V_sample_VUF(nc+idxB, :);
    Vc_gmm = V_sample_VUF(idxC, :) + 1i * V_sample_VUF(nc+idxC, :);

    % 序分量变换
    V_seq_gmm = T_seq * [Va_gmm; Vb_gmm; Vc_gmm];  % [3 x nSamples]

    Vmag_pos_gmm = abs(V_seq_gmm(2, :));  % 正序幅值
    Vmag_neg_gmm = abs(V_seq_gmm(3, :));  % 负序幅值

    VUF_gmm_samples = Vmag_neg_gmm ./ Vmag_pos_gmm * 100;

    mean_VUF_exact(b) = mean(VUF_gmm_samples);
    var_VUF_exact(b)  = var(VUF_gmm_samples);
end
clear V_sample_VUF;  % 释放内存
fprintf('[VUF] 精确矩计算完成\n');

% 10.8 近似计算VUF均值方差（改进的二次型比值Delta法）
%     改进要点:
%       (1) 直接在|V-|²/|V+|²比值上做Delta法，避免两步近似累积误差
%       (2) 纳入Cov(|V+|², |V-|²)项（二者共享电压向量，存在相关性）
%       (3) 再对√W做二阶校正: E[√W] ≈ √E[W]·(1 - Var(W)/(8·E[W]²))
fprintf('[VUF] 改进的比值Delta法计算近似矩...\n');

mean_VUF_approx = zeros(n3ph, 1);
var_VUF_approx  = zeros(n3ph, 1);

for b = 1:n3ph
    % Step 1: 各GMM分量的矩（从GNCS参数解析计算）
    mQ_neg = zeros(gmm_num, 1);  vQ_neg = zeros(gmm_num, 1);
    mQ_pos = zeros(gmm_num, 1);  vQ_pos = zeros(gmm_num, 1);
    covQ   = zeros(gmm_num, 1);

    for k = 1:gmm_num
        % 负序 |V-|² 的均值和方差
        mQ_neg(k) = Vneg2_const{b,k} + sum(Vneg2_lumda{b,k}.*(1+Vneg2_delta{b,k}));
        vQ_neg(k) = 2*sum(Vneg2_lumda{b,k}.^2.*(1+2*Vneg2_delta{b,k}));

        % 正序 |V+|² 的均值和方差
        mQ_pos(k) = Vpos2_const{b,k} + sum(Vpos2_lumda{b,k}.*(1+Vpos2_delta{b,k}));
        vQ_pos(k) = 2*sum(Vpos2_lumda{b,k}.^2.*(1+2*Vpos2_delta{b,k}));

        % Cov_k(|V+|², |V-|²) = 2·tr(M+·Σ_k·M-·Σ_k) + 4·μ_k'·M+·Σ_k·M-·μ_k
        mu_k = gmm_fit.mu(k,:)';
        Sigma_k = gmm_fit.Sigma(:,:,k);
        MSp = M_Vpos2{b} * Sigma_k;
        MSn = M_Vneg2{b} * Sigma_k;
        covQ(k) = 2*trace(MSp * MSn) + 4*mu_k'*M_Vpos2{b}*Sigma_k*M_Vneg2{b}*mu_k;
    end

    % Step 2: GMM混合矩
    E_Qneg = mQ_neg' * com_pi';
    E_Qpos = mQ_pos' * com_pi';
    V_Qneg = vQ_neg' * com_pi' + (mQ_neg.^2)' * com_pi' - E_Qneg^2;
    V_Qpos = vQ_pos' * com_pi' + (mQ_pos.^2)' * com_pi' - E_Qpos^2;
    Cov_QnQp = covQ' * com_pi' + (mQ_neg.*mQ_pos)' * com_pi' - E_Qneg*E_Qpos;

    % Step 3: 比值 W = |V-|²/|V+|² 的Delta法
    %   E[W] ≈ E[Qn]/E[Qp] · (1 + V(Qp)/E[Qp]² - Cov/(E[Qn]·E[Qp]))
    %   Var(W) ≈ (E[Qn]/E[Qp])² · (V(Qn)/E[Qn]² + V(Qp)/E[Qp]² - 2·Cov/(E[Qn]·E[Qp]))
    r0 = E_Qneg / E_Qpos;
    cv2_n = V_Qneg / E_Qneg^2;
    cv2_p = V_Qpos / E_Qpos^2;
    rho_np = Cov_QnQp / (E_Qneg * E_Qpos);

    E_W = r0 * (1 + cv2_p - rho_np);
    V_W = r0^2 * (cv2_n + cv2_p - 2*rho_np);

    % Step 4: VUF = 100·√W
    %   E[√W] ≈ √(E[W])·(1 - V(W)/(8·E[W]²))
    %   Var(√W) ≈ V(W)/(4·E[W])
    cv2_W = V_W / E_W^2;
    mean_VUF_approx(b) = 100 * sqrt(E_W) * (1 - cv2_W/8);
    var_VUF_approx(b)  = 100^2 * V_W / (4 * E_W);
end
fprintf('[VUF] 近似矩计算完成\n');

% 10.9 误差对比
error_VUF_mean_exact  = abs((mean_VUF_mc - mean_VUF_exact)  ./ mean_VUF_mc) * 100;
error_VUF_mean_approx = abs((mean_VUF_mc - mean_VUF_approx) ./ mean_VUF_mc) * 100;
error_VUF_var_exact   = abs((var_VUF_mc  - var_VUF_exact)   ./ var_VUF_mc)  * 100;
error_VUF_var_approx  = abs((var_VUF_mc  - var_VUF_approx)  ./ var_VUF_mc)  * 100;

fprintf('[VUF] 均值精确最大误差: %.4f%%, 近似最大误差: %.4f%%\n', ...
    max(error_VUF_mean_exact), max(error_VUF_mean_approx));
fprintf('[VUF] 方差精确最大误差: %.4f%%, 近似最大误差: %.4f%%\n', ...
    max(error_VUF_var_exact), max(error_VUF_var_approx));

% 打印各节点详细VUF误差
fprintf('\n  %-12s  %10s  %10s  %10s  %10s\n', ...
    '三相母线', '均值精确%', '均值近似%', '方差精确%', '方差近似%');
fprintf('  %s\n', repmat('-', 1, 56));
for b = 1:n3ph
    fprintf('  %-12s  %10.4f  %10.4f  %10.4f  %10.4f\n', ...
        threePhaseBuses{b}, error_VUF_mean_exact(b), error_VUF_mean_approx(b), ...
        error_VUF_var_exact(b), error_VUF_var_approx(b));
end

%% ========================================================================
%  段落11: 序电压分析 (Fig.11-12) 与 VUF图 (Fig.13)
%  ========================================================================
fprintf('\n[序电压] 计算序电压精确矩 (GNCS采样)...\n');

% 11.1 序电压精确矩 (GNCS采样)
n_seq_samp = 1e5;
mean_Vpos_exact = zeros(n3ph,1); var_Vpos_exact = zeros(n3ph,1);
mean_Vneg_exact = zeros(n3ph,1); var_Vneg_exact = zeros(n3ph,1);
mean_V0_exact   = zeros(n3ph,1); var_V0_exact   = zeros(n3ph,1);

for b = 1:n3ph
    % 正序 |V+|
    smp = [];
    for k = 1:gmm_num
        nk = round(n_seq_samp*com_pi(k)); if nk==0, continue; end
        Z = randn(length(Vpos2_lumda{b,k}), nk);
        sq = Vpos2_const{b,k} + sum(Vpos2_lumda{b,k}.*(Z+sqrt(Vpos2_delta{b,k})).^2,1);
        smp = [smp, sqrt(max(sq,0))];
    end
    mean_Vpos_exact(b) = mean(smp); var_Vpos_exact(b) = var(smp);

    % 负序 |V-|
    smp = [];
    for k = 1:gmm_num
        nk = round(n_seq_samp*com_pi(k)); if nk==0, continue; end
        Z = randn(length(Vneg2_lumda{b,k}), nk);
        sq = Vneg2_const{b,k} + sum(Vneg2_lumda{b,k}.*(Z+sqrt(Vneg2_delta{b,k})).^2,1);
        smp = [smp, sqrt(max(sq,0))];
    end
    mean_Vneg_exact(b) = mean(smp); var_Vneg_exact(b) = var(smp);

    % 零序 |V0|
    smp = [];
    for k = 1:gmm_num
        nk = round(n_seq_samp*com_pi(k)); if nk==0, continue; end
        Z = randn(length(V0sq_lumda{b,k}), nk);
        sq = V0sq_const{b,k} + sum(V0sq_lumda{b,k}.*(Z+sqrt(V0sq_delta{b,k})).^2,1);
        smp = [smp, sqrt(max(sq,0))];
    end
    mean_V0_exact(b) = mean(smp); var_V0_exact(b) = var(smp);
end

% 序电压误差（使用中位数阈值避免小分母效应）
thr_mean_Vpos = max(abs(mean_Vpos_mc)) * 0.01;
thr_mean_Vneg = max(abs(mean_Vneg_mc)) * 0.01;
thr_mean_V0   = max(abs(mean_V0_mc))   * 0.01;
thr_var_Vpos  = max(abs(var_Vpos_mc))  * 0.01;
thr_var_Vneg  = max(abs(var_Vneg_mc))  * 0.01;
thr_var_V0    = max(abs(var_V0_mc))    * 0.01;

err_Vpos_mean = abs((mean_Vpos_mc - mean_Vpos_exact) ./ max(abs(mean_Vpos_mc), thr_mean_Vpos))*100;
err_Vneg_mean = abs((mean_Vneg_mc - mean_Vneg_exact) ./ max(abs(mean_Vneg_mc), thr_mean_Vneg))*100;
err_V0_mean   = abs((mean_V0_mc - mean_V0_exact)     ./ max(abs(mean_V0_mc),   thr_mean_V0))*100;
err_Vpos_var  = abs((var_Vpos_mc - var_Vpos_exact)    ./ max(abs(var_Vpos_mc),  thr_var_Vpos))*100;
err_Vneg_var  = abs((var_Vneg_mc - var_Vneg_exact)    ./ max(abs(var_Vneg_mc),  thr_var_Vneg))*100;
err_V0_var    = abs((var_V0_mc - var_V0_exact)        ./ max(abs(var_V0_mc),    thr_var_V0))*100;

fprintf('[序电压] 正序均值最大误差: %.4f%%, 负序: %.4f%%, 零序: %.4f%%\n', ...
    max(err_Vpos_mean), max(err_Vneg_mean), max(err_V0_mean));

% 三相节点标签（按最小节点编号从小到大排序，使用母线名）
[~, sortIdx_3ph] = sort(min(threePhaseBusIdx, [], 2));
busLabels_3ph = cell(n3ph, 1);
for b = 1:n3ph
    sb = sortIdx_3ph(b);
    busLabels_3ph{b} = busLabels_of_node{threePhaseBusIdx(sb,1)};
end
busLabels_3ph_short = busLabels_3ph;
% 对误差向量也按相同顺序排列
err_Vpos_mean = err_Vpos_mean(sortIdx_3ph);
err_Vneg_mean = err_Vneg_mean(sortIdx_3ph);
err_V0_mean   = err_V0_mean(sortIdx_3ph);
err_Vpos_var  = err_Vpos_var(sortIdx_3ph);
err_Vneg_var  = err_Vneg_var(sortIdx_3ph);
err_V0_var    = err_V0_var(sortIdx_3ph);

% (旧 Fig.12 序电压PDF已删除，移至段落12.7新图5)


%% ========================================================================
%  段落12: SOTA方法对比 (PEM / CM)
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  段落12: SOTA方法对比 (PEM / CM)\n');
fprintf('%s\n', repmat('=',1,70));

% --- 12.0 缓存检查 ---
compCacheFile = fullfile(cacheDir, sprintf('Comparison_%s_v2.mat', caseName));

if isfile(compCacheFile)
    fprintf('[对比] 加载缓存: %s\n', compCacheFile);
    compData = load(compCacheFile);
    V_pem_all     = compData.V_pem_all;
    w_pem         = compData.w_pem;
    w0_pem        = compData.w0_pem;
    V_pem_center  = compData.V_pem_center;
    V_cum_base    = compData.V_cum_base;
    J_v           = compData.J_v;
    t_pem_dss     = compData.t_pem_dss;
    t_cum_dss     = compData.t_cum_dss;
    mu_X          = compData.mu_X;
    sigma_X       = compData.sigma_X;
    X_data        = compData.X_data;
    m_rv          = compData.m_rv;
    fprintf('[对比] 缓存加载完成\n');
else
    % --- 12.1 OpenDSS初始化 ---
    fprintf('[对比] 初始化 OpenDSS COM...\n');
    DSSObj = actxserver('OpenDSSEngine.DSS');
    assert(DSSObj.Start(0), '[OpenDSS] 启动失败');

    % 输入随机变量统计量
    m_rv = 24;  % 负荷缩放曲线数
    X_data = zeros(N, m_rv);
    for c = 1:m_rv
        X_data(:, c) = load_data{:, c};
    end
    mu_X    = mean(X_data);      % [1 x m]
    sigma_X = std(X_data);       % [1 x m]
    skew_X  = skewness(X_data);  % [1 x m]

    fprintf('[对比] 输入维度 m=%d, 样本数 N=%d\n', m_rv, N);

    % =====================================================================
    % 12.2 PEM (Hong's 2m+1)
    % =====================================================================
    fprintf('\n[PEM] Hong 2m+1 点估计法 (共 %d 次潮流)...\n', 2*m_rv+1);

    % 计算集中点和权重
    xi1 = zeros(m_rv, 1); xi2 = zeros(m_rv, 1);
    for i = 1:m_rv
        lam3 = skew_X(i);
        xi1(i) = lam3/2 + sqrt(m_rv + (lam3/2)^2);
        xi2(i) = lam3/2 - sqrt(m_rv + (lam3/2)^2);
    end

    w1 = zeros(m_rv, 1); w2 = zeros(m_rv, 1);
    for i = 1:m_rv
        w1(i) = (-xi2(i)) / (m_rv * xi1(i) * (xi1(i) - xi2(i)));
        w2(i) =   xi1(i)  / (m_rv * xi2(i) * (xi2(i) - xi1(i)));
    end
    w0_pem = 1 - sum(w1 + w2);

    % 运行中心点
    tic;
    V_pem_center = tools.eval_opendss_single(DSSObj, dssFile, mu_X(:), ...
        load_sequence(:), loadKW_base, loadKvar_base, nLoad);

    % 运行2m集中点
    V_pem_all = zeros(nc, 2*m_rv);
    w_pem = zeros(2*m_rv, 1);
    h_pem = waitbar(0, 'PEM: 0/48', 'Name', 'PEM 2m+1');
    for i = 1:m_rv
        x_pt1 = mu_X(:); x_pt1(i) = mu_X(i) + xi1(i)*sigma_X(i);
        x_pt2 = mu_X(:); x_pt2(i) = mu_X(i) + xi2(i)*sigma_X(i);

        V_pem_all(:, 2*i-1) = tools.eval_opendss_single(DSSObj, dssFile, ...
            x_pt1, load_sequence(:), loadKW_base, loadKvar_base, nLoad);
        V_pem_all(:, 2*i)   = tools.eval_opendss_single(DSSObj, dssFile, ...
            x_pt2, load_sequence(:), loadKW_base, loadKvar_base, nLoad);
        w_pem(2*i-1) = w1(i);
        w_pem(2*i)   = w2(i);
        waitbar(i/m_rv, h_pem, sprintf('PEM: %d/%d', 2*i, 2*m_rv));
    end
    close(h_pem);
    t_pem_dss = toc;
    fprintf('[PEM] OpenDSS评估完成, 耗时: %.2f s\n', t_pem_dss);

    % =====================================================================
    % 12.3 CM (半不变量法 / 1阶线性化)
    % =====================================================================
    fprintf('\n[CM] 数值Jacobian计算 (共 %d 次潮流)...\n', 2*m_rv+1);

    % 基准点（均值点）— 复用PEM中心点
    V_cum_base = V_pem_center;

    % 数值Jacobian: ∂V_complex/∂x_i (中心差分)
    delta_frac = 0.01;  % 扰动量 = 1% 标准差
    J_v = zeros(2*nc, m_rv);  % 实同构Jacobian [2*nc x m]
    tic;
    h_cum = waitbar(0, 'CM: 0/48', 'Name', 'CM Jacobian');
    for i = 1:m_rv
        delta_i = delta_frac * sigma_X(i);
        x_plus  = mu_X(:); x_plus(i)  = mu_X(i) + delta_i;
        x_minus = mu_X(:); x_minus(i) = mu_X(i) - delta_i;

        V_plus  = tools.eval_opendss_single(DSSObj, dssFile, x_plus, ...
            load_sequence(:), loadKW_base, loadKvar_base, nLoad);
        V_minus = tools.eval_opendss_single(DSSObj, dssFile, x_minus, ...
            load_sequence(:), loadKW_base, loadKvar_base, nLoad);

        % 实同构形式的Jacobian
        dV_re = (real(V_plus) - real(V_minus)) / (2*delta_i);
        dV_im = (imag(V_plus) - imag(V_minus)) / (2*delta_i);
        J_v(:, i) = [dV_re; dV_im];

        waitbar(i/m_rv, h_cum, sprintf('CM: %d/%d', 2*i, 2*m_rv));
    end
    close(h_cum);
    t_cum_dss = toc;
    fprintf('[CM] Jacobian计算完成, 耗时: %.2f s\n', t_cum_dss);

    % --- 保存缓存 ---
    fprintf('[对比] 保存缓存: %s\n', compCacheFile);
    save(compCacheFile, 'V_pem_all', 'w_pem', 'w0_pem', 'V_pem_center', ...
        'V_cum_base', 'J_v', ...
        't_pem_dss', 't_cum_dss', ...
        'mu_X', 'sigma_X', 'X_data', 'm_rv', '-v7.3');
    fprintf('[对比] 缓存保存完成\n');
    DSSObj.delete;
    clear DSSObj;
    fprintf('[对比] OpenDSS COM已释放\n');
end

% =====================================================================
% 12.5 后处理: 计算各方法的线路功率和VUF
% =====================================================================
fprintf('\n[对比] 后处理: 计算线路功率和VUF...\n');
tic;

Sigma_X = cov(X_data);  % [m x m] 输入协方差矩阵

% --- (A) PEM: 加权矩计算 ---
nSP = length(samePhaseIdx);
v_center = [real(V_pem_center); imag(V_pem_center)];  % [2*nc x 1]

% 线路功率
mean_P_pem = zeros(nBranch, 1); var_P_pem = zeros(nBranch, 1);
mean_Q_pem = zeros(nBranch, 1); var_Q_pem = zeros(nBranch, 1);
% PEM第3矩 (用于Gram-Charlier)
skew_P_pem = zeros(nBranch, 1); skew_Q_pem = zeros(nBranch, 1);

for idx = 1:nSP
    l = samePhaseIdx(idx);
    P_c = v_center' * Line_P_matrix{l} * v_center;
    Q_c = v_center' * Line_Q_matrix{l} * v_center;

    P_pts = zeros(2*m_rv, 1); Q_pts = zeros(2*m_rv, 1);
    for k = 1:2*m_rv
        vk = [real(V_pem_all(:,k)); imag(V_pem_all(:,k))];
        P_pts(k) = vk' * Line_P_matrix{l} * vk;
        Q_pts(k) = vk' * Line_Q_matrix{l} * vk;
    end

    mean_P_pem(l) = w0_pem*P_c + sum(w_pem .* P_pts);
    mean_Q_pem(l) = w0_pem*Q_c + sum(w_pem .* Q_pts);
    EP2 = w0_pem*P_c^2 + sum(w_pem .* P_pts.^2);
    EQ2 = w0_pem*Q_c^2 + sum(w_pem .* Q_pts.^2);
    var_P_pem(l) = max(EP2 - mean_P_pem(l)^2, 0);
    var_Q_pem(l) = max(EQ2 - mean_Q_pem(l)^2, 0);
    EP3 = w0_pem*P_c^3 + sum(w_pem .* P_pts.^3);
    EQ3 = w0_pem*Q_c^3 + sum(w_pem .* Q_pts.^3);
    mu3_P = EP3 - 3*mean_P_pem(l)*EP2 + 2*mean_P_pem(l)^3;
    mu3_Q = EQ3 - 3*mean_Q_pem(l)*EQ2 + 2*mean_Q_pem(l)^3;
    skew_P_pem(l) = mu3_P / max(var_P_pem(l)^1.5, eps);
    skew_Q_pem(l) = mu3_Q / max(var_Q_pem(l)^1.5, eps);
end

% VUF
mean_VUF_pem = zeros(n3ph, 1); var_VUF_pem = zeros(n3ph, 1);
skew_VUF_pem = zeros(n3ph, 1);
for b = 1:n3ph
    idxA = threePhaseBusIdx(b,1); idxB = threePhaseBusIdx(b,2); idxC = threePhaseBusIdx(b,3);

    Vc = [V_pem_center(idxA); V_pem_center(idxB); V_pem_center(idxC)];
    Vs = T_seq * Vc;
    vuf_c = abs(Vs(3))/abs(Vs(2))*100;

    vuf_pts = zeros(2*m_rv, 1);
    for k = 1:2*m_rv
        Vk = [V_pem_all(idxA,k); V_pem_all(idxB,k); V_pem_all(idxC,k)];
        Vsk = T_seq * Vk;
        vuf_pts(k) = abs(Vsk(3))/abs(Vsk(2))*100;
    end

    mean_VUF_pem(b) = w0_pem*vuf_c + sum(w_pem.*vuf_pts);
    EV2 = w0_pem*vuf_c^2 + sum(w_pem.*vuf_pts.^2);
    var_VUF_pem(b) = max(EV2 - mean_VUF_pem(b)^2, 0);
    EV3 = w0_pem*vuf_c^3 + sum(w_pem.*vuf_pts.^3);
    mu3_v = EV3 - 3*mean_VUF_pem(b)*EV2 + 2*mean_VUF_pem(b)^3;
    skew_VUF_pem(b) = mu3_v / max(var_VUF_pem(b)^1.5, eps);
end
% PEM: 各节点电压幅值 |V| 的加权均值/方差/偏度 (用于Fig.4a/4b及PDF对比)
mean_Vmag_pem = zeros(nc, 1); var_Vmag_pem = zeros(nc, 1); skew_Vmag_pem = zeros(nc, 1);
for i = 1:nc
    vmag_c = abs(V_pem_center(i));
    vmag_pts = abs(V_pem_all(i, :)).';  % [2*m_rv x 1]
    mean_Vmag_pem(i) = w0_pem * vmag_c + sum(w_pem .* vmag_pts);
    EV2 = w0_pem * vmag_c^2 + sum(w_pem .* vmag_pts.^2);
    var_Vmag_pem(i) = max(EV2 - mean_Vmag_pem(i)^2, 0);
    EV3 = w0_pem * vmag_c^3 + sum(w_pem .* vmag_pts.^3);
    mu3 = EV3 - 3*mean_Vmag_pem(i)*EV2 + 2*mean_Vmag_pem(i)^3;
    skew_Vmag_pem(i) = mu3 / max(var_Vmag_pem(i)^1.5, eps);
end
fprintf('[PEM] 后处理完成\n');

% --- (B) CM (Cumulant): 线性化传播 ---
v_base = [real(V_cum_base); imag(V_cum_base)];  % [2*nc x 1]

% 将MC输入通过线性化模型传播: v ≈ v_base + J_v · (x - μ)
dX = X_data - mu_X;  % [N x m]
V_cum_all = v_base + J_v * dX';  % [2*nc x N] 线性化电压

mean_P_cum = zeros(nBranch, 1); var_P_cum = zeros(nBranch, 1);
mean_Q_cum = zeros(nBranch, 1); var_Q_cum = zeros(nBranch, 1);
skew_P_cum = zeros(nBranch, 1); skew_Q_cum = zeros(nBranch, 1);

for idx = 1:nSP
    l = samePhaseIdx(idx);
    P_cum_l = dot(V_cum_all, Line_P_matrix{l} * V_cum_all)';  % [N x 1]
    Q_cum_l = dot(V_cum_all, Line_Q_matrix{l} * V_cum_all)';
    mean_P_cum(l) = mean(P_cum_l);
    var_P_cum(l)  = var(P_cum_l);
    mean_Q_cum(l) = mean(Q_cum_l);
    var_Q_cum(l)  = var(Q_cum_l);
    skew_P_cum(l) = skewness(P_cum_l);
    skew_Q_cum(l) = skewness(Q_cum_l);
end

% VUF (从线性化电压计算)
V_cum_complex = V_cum_all(1:nc, :) + 1i * V_cum_all(nc+1:end, :);  % [nc x N]
mean_VUF_cum = zeros(n3ph, 1); var_VUF_cum = zeros(n3ph, 1);
skew_VUF_cum = zeros(n3ph, 1);
for b = 1:n3ph
    idxA = threePhaseBusIdx(b,1); idxB = threePhaseBusIdx(b,2); idxC = threePhaseBusIdx(b,3);
    V_abc_cum = [V_cum_complex(idxA,:); V_cum_complex(idxB,:); V_cum_complex(idxC,:)];
    V_seq_cum = T_seq * V_abc_cum;
    VUF_cum = abs(V_seq_cum(3,:)) ./ abs(V_seq_cum(2,:)) * 100;
    mean_VUF_cum(b) = mean(VUF_cum);
    var_VUF_cum(b)  = var(VUF_cum);
    skew_VUF_cum(b) = skewness(VUF_cum(:));
end
% CM: 各节点电压幅值 |V| 的均值/方差/偏度
Vmag_cum = abs(V_cum_complex);  % [nc x N]
mean_Vmag_cum = mean(Vmag_cum, 2);
var_Vmag_cum  = var(Vmag_cum, 0, 2);
skew_Vmag_cum = zeros(nc, 1);
for i = 1:nc
    skew_Vmag_cum(i) = skewness(Vmag_cum(i,:)');
end
fprintf('[CM] 后处理完成\n');

t_post = toc;
fprintf('[对比] 后处理总耗时: %.2f s\n', t_post);

% =====================================================================
% 12.6 误差计算
% =====================================================================
fprintf('\n[对比] 计算各方法误差...\n');

% MC基准 (同相线路)
mean_P_mc_sp = mean_Line_P_mc(samePhaseIdx);
mean_Q_mc_sp = mean_Line_Q_mc(samePhaseIdx);
var_P_mc_sp  = var_Line_P_mc(samePhaseIdx);
var_Q_mc_sp  = var_Line_Q_mc(samePhaseIdx);

% 线路功率均值误差 (%)
thr_P = 0.01 * median(abs(mean_P_mc_sp));
thr_Q = 0.01 * median(abs(mean_Q_mc_sp));
err_mean_P = @(x) abs(mean_P_mc_sp - x(samePhaseIdx)) ./ max(abs(mean_P_mc_sp), thr_P) * 100;
err_mean_Q = @(x) abs(mean_Q_mc_sp - x(samePhaseIdx)) ./ max(abs(mean_Q_mc_sp), thr_Q) * 100;

errMP_aqppf = err_mean_P(mean_Line_P_this);
errMP_pem   = err_mean_P(mean_P_pem);
errMP_cum   = err_mean_P(mean_P_cum);

errMQ_aqppf = err_mean_Q(mean_Line_Q_this);
errMQ_pem   = err_mean_Q(mean_Q_pem);
errMQ_cum   = err_mean_Q(mean_Q_cum);

% 线路功率方差误差 (%)
thr_vP = 0.01 * median(abs(var_P_mc_sp));
thr_vQ = 0.01 * median(abs(var_Q_mc_sp));
err_var_P = @(x) abs(var_P_mc_sp - x(samePhaseIdx)) ./ max(abs(var_P_mc_sp), thr_vP) * 100;
err_var_Q = @(x) abs(var_Q_mc_sp - x(samePhaseIdx)) ./ max(abs(var_Q_mc_sp), thr_vQ) * 100;

errVP_aqppf = err_var_P(var_Line_P_this);
errVP_pem   = err_var_P(var_P_pem);
errVP_cum   = err_var_P(var_P_cum);

errVQ_aqppf = err_var_Q(var_Line_Q_this);
errVQ_pem   = err_var_Q(var_Q_pem);
errVQ_cum   = err_var_Q(var_Q_cum);

% VUF误差 (%)
errVUF_mean_aqppf = abs((mean_VUF_mc - mean_VUF_exact)  ./ mean_VUF_mc) * 100;
errVUF_mean_pem   = abs((mean_VUF_mc - mean_VUF_pem)    ./ mean_VUF_mc) * 100;
errVUF_mean_cum   = abs((mean_VUF_mc - mean_VUF_cum)    ./ mean_VUF_mc) * 100;

errVUF_var_aqppf  = abs((var_VUF_mc - var_VUF_exact)    ./ var_VUF_mc) * 100;
errVUF_var_pem    = abs((var_VUF_mc - var_VUF_pem)      ./ var_VUF_mc) * 100;
errVUF_var_cum    = abs((var_VUF_mc - var_VUF_cum)      ./ var_VUF_mc) * 100;

% 电压幅值误差 (%)
thr_vmag_m = max(mean_Vmag_mc) * 1e-6 + 1e-12;
thr_vmag_v = max(var_Vmag_mc) * 1e-6 + 1e-12;
errVmag_mean_aqppf = error_Vmag_mean_exact;
errVmag_mean_pem   = abs((mean_Vmag_mc - mean_Vmag_pem)   ./ max(mean_Vmag_mc, thr_vmag_m)) * 100;
errVmag_mean_cum   = abs((mean_Vmag_mc - mean_Vmag_cum)   ./ max(mean_Vmag_mc, thr_vmag_m)) * 100;
errVmag_var_aqppf  = error_Vmag_var_exact;
errVmag_var_pem    = abs((var_Vmag_mc - var_Vmag_pem)    ./ max(var_Vmag_mc, thr_vmag_v)) * 100;
errVmag_var_cum    = abs((var_Vmag_mc - var_Vmag_cum)    ./ max(var_Vmag_mc, thr_vmag_v)) * 100;

fprintf('[对比] 线路P均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errMP_aqppf), max(errMP_pem), max(errMP_cum));
fprintf('[对比] 线路P方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVP_aqppf), max(errVP_pem), max(errVP_cum));
fprintf('[对比] VUF均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVUF_mean_aqppf), max(errVUF_mean_pem), max(errVUF_mean_cum));
fprintf('[对比] VUF方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVUF_var_aqppf), max(errVUF_var_pem), max(errVUF_var_cum));
fprintf('[对比] 电压幅值均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVmag_mean_aqppf), max(errVmag_mean_pem), max(errVmag_mean_cum));
fprintf('[对比] 电压幅值方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVmag_var_aqppf), max(errVmag_var_pem), max(errVmag_var_cum));

% =====================================================================
% 12.7 对比图表 (目标产物: Fig3, Fig4, Tab1, Fig5, Fig6, Tab2)
% =====================================================================
fprintf('\n[对比] 生成对比图表...\n');

method_names = {'本文', 'PEM', 'CM'};
method_colors = {[0 0.45 0.74], [0.85 0.33 0.1], [0.49 0.18 0.56]};
fontFig = 22;

% === Fig.3: 代表性节点电压幅值PDF对比 (MCS/本文/PEM/CM) ===
fprintf('\n[图3] 电压幅值PDF对比...\n');
show_vmag_labels = {'611-C', '646-B'};
show_nodes_vmag = zeros(1, length(show_vmag_labels));
for idx = 1:length(show_vmag_labels)
    match = find(strcmp(nodeLabels, show_vmag_labels{idx}));
    assert(~isempty(match), '未找到节点 %s', show_vmag_labels{idx});
    show_nodes_vmag(idx) = match(1);
end
for idx = 1:length(show_nodes_vmag)
    ni = show_nodes_vmag(idx);
    Vb = V_base_node(ni);
    Vmag_mc_pu = Vmag_mc(:, ni) / Vb;
    x_min = min(Vmag_mc_pu) * 0.98; x_max = max(Vmag_mc_pu) * 1.02;

    figure('Units','pixels','Position',[100+60*idx, 100+60*idx, 700, 500]);
    histogram(Vmag_mc_pu, 80, 'Normalization','pdf', 'FaceAlpha',0.3, ...
        'FaceColor',[0.07 0.44 0.75], 'EdgeColor','none');
    hold on;
    xr_pu = linspace(x_min, x_max, 300);

    % 本文: Davies解析PDF  f_{|V|}(v) = 2v·f_W(v^2)
    pv_aqppf = zeros(size(xr_pu));
    for k = 1:gmm_num
        if isempty(Vmag2_lumda{ni,k}), continue; end
        w_vals = xr_pu.^2 * Vb^2;
        f_W = gx2pdf_davies(w_vals, Vmag2_lumda{ni,k}', ...
            ones(1,length(Vmag2_lumda{ni,k})), Vmag2_delta{ni,k}', ...
            Vmag2_const{ni,k}, 0);
        f_pu_k = 2 * xr_pu * Vb^2 .* f_W;
        pv_aqppf = pv_aqppf + com_pi(k) * f_pu_k;
    end
    plot(xr_pu, pv_aqppf, '-', 'Color', method_colors{1}, 'LineWidth', 4);

    % PEM: Gram-Charlier
    mu_gc = mean_Vmag_pem(ni)/Vb; sig_gc = sqrt(max(var_Vmag_pem(ni),0))/Vb;
    z_gc = (xr_pu - mu_gc) / max(sig_gc, eps);
    H3 = z_gc.^3 - 3*z_gc;
    f_pem = normpdf(z_gc)/max(sig_gc,eps) .* (1 + skew_Vmag_pem(ni)/6*H3);
    f_pem = max(f_pem, 0);
    plot(xr_pu, f_pem, '--', 'Color', method_colors{2}, 'LineWidth', 4);

    % CM: Gram-Charlier
    mu_gc = mean_Vmag_cum(ni)/Vb; sig_gc = sqrt(max(var_Vmag_cum(ni),0))/Vb;
    z_gc = (xr_pu - mu_gc) / max(sig_gc, eps);
    H3 = z_gc.^3 - 3*z_gc;
    f_cum = normpdf(z_gc)/max(sig_gc,eps) .* (1 + skew_Vmag_cum(ni)/6*H3);
    f_cum = max(f_cum, 0);
    plot(xr_pu, f_cum, ':', 'Color', method_colors{3}, 'LineWidth', 4);

    hold off;
    set(gca,'fontsize',fontFig,'fontname',fontName_en,'box','off');
    if strcmp(show_vmag_labels{idx}, '611-C')
        xlim([0.96, 1.02]);
    elseif strcmp(show_vmag_labels{idx}, '646-B')
        xlim([0.985, 1.015]);
    end
    xlabel('电压幅值 (p.u.)','FontName',fontName_cn,'FontSize',fontFig);
    ylabel('PDF','FontSize',fontFig);
    legend({'MCS','本文','PEM','CM'},'FontName',fontName_cn,'FontSize',fontFig,'Location','best');
    title(sprintf('%s', nodeLabels{ni}),'FontName',fontName_cn,'FontSize',fontFig);
    saveas(gcf, fullfile(resultDir, sprintf('Fig03_%s_VmagPDF.png', nodeLabels{ni})));
end
fprintf('[图3] 已保存\n');

% === Fig.4: 支路 645→632 (两相B+C) 功率PDF对比，4张独立图 ===
fprintf('\n[图4] 支路645→632功率PDF对比（4张独立图）...\n');

% 线路键按字典序: '632-645'（busFrom='645'>busTo='632'，故key=busTo-busFrom）
target_line_key = '632-645';
if ~physLineMap.isKey(target_line_key)
    warning('[图4] 未找到支路 %s，跳过。', target_line_key);
else
    info_t = physLineMap(target_line_key);
    act_ph_t = find(info_t.brIdx > 0);  % 应为 [2,3] 即B、C相

    % 只绘制有功功率 P，B相和C相各一张
    subplot_specs = {};
    for phi_idx = 1:length(act_ph_t)
        subplot_specs{end+1} = {act_ph_t(phi_idx), 'P'};
    end

    fig4_suffixes = {'B_P', 'C_P'};

    for si = 1:length(subplot_specs)
        ph    = subplot_specs{si}{1};
        pq    = subplot_specs{si}{2};
        l     = info_t.brIdx(ph);

        figure('Units','pixels','Position',[100+80*si, 100+80*si, 560, 420]);
        hold on;

        if strcmp(pq, 'P')
            % --- 有功功率 ---
            P_mc_MW = Line_Power_P_MC{l} / 1e6;
            histogram(P_mc_MW, 60, 'Normalization','pdf', 'FaceAlpha',0.3, ...
                'FaceColor',[0.07 0.44 0.75], 'EdgeColor','none');

            xr = linspace(min(Line_Power_P_MC{l}), max(Line_Power_P_MC{l}), n_pdf_pts);
            pv = zeros(1, n_pdf_pts);
            for k = 1:gmm_num
                pv = pv + com_pi(k)*gx2pdf_davies(xr, ...
                    Line_P_lumda{l,k}', ones(1,length(Line_P_lumda{l,k})), ...
                    Line_P_delta{l,k}', Line_P_const{l,k}, 0);
            end
            plot(xr/1e6, pv*1e6, '-', 'Color', method_colors{1}, 'LineWidth', 2.5);

            mu_p = mean_P_pem(l); sig_p = sqrt(max(var_P_pem(l),0)); gam_p = skew_P_pem(l);
            z_p = (xr - mu_p) / max(sig_p, eps);
            f_pem = normpdf(z_p)/max(sig_p,eps) .* (1 + gam_p/6*(z_p.^3 - 3*z_p));
            plot(xr/1e6, max(f_pem,0)*1e6, '--', 'Color', method_colors{2}, 'LineWidth', 2);

            mu_cp = mean_P_cum(l); sig_cp = sqrt(max(var_P_cum(l),0)); gam_cp = skew_P_cum(l);
            z_cp = (xr - mu_cp) / max(sig_cp, eps);
            f_cum = normpdf(z_cp)/max(sig_cp,eps) .* (1 + gam_cp/6*(z_cp.^3 - 3*z_cp));
            plot(xr/1e6, max(f_cum,0)*1e6, ':', 'Color', method_colors{3}, 'LineWidth', 2);

            xlabel('有功功率 (MW)', 'FontName',fontName_cn, 'FontSize',fontSZ);
            title(sprintf('%s\x2192%s %s相 P', info_t.from_bus, info_t.to_bus, phaseChar_arr{ph}), ...
                'FontName',fontName_cn, 'FontSize',fontSZ);
        else
            % --- 无功功率 ---
            Q_mc_MV = Line_Power_Q_MC{l} / 1e6;
            histogram(Q_mc_MV, 60, 'Normalization','pdf', 'FaceAlpha',0.3, ...
                'FaceColor',[0.07 0.44 0.75], 'EdgeColor','none');

            xr = linspace(min(Line_Power_Q_MC{l}), max(Line_Power_Q_MC{l}), n_pdf_pts);
            pv = zeros(1, n_pdf_pts);
            for k = 1:gmm_num
                pv = pv + com_pi(k)*gx2pdf_davies(xr, ...
                    Line_Q_lumda{l,k}', ones(1,length(Line_Q_lumda{l,k})), ...
                    Line_Q_delta{l,k}', Line_Q_const{l,k}, 0);
            end
            plot(xr/1e6, pv*1e6, '-', 'Color', method_colors{1}, 'LineWidth', 2.5);

            mu_q = mean_Q_pem(l); sig_q = sqrt(max(var_Q_pem(l),0)); gam_q = skew_Q_pem(l);
            z_q = (xr - mu_q) / max(sig_q, eps);
            f_pem = normpdf(z_q)/max(sig_q,eps) .* (1 + gam_q/6*(z_q.^3 - 3*z_q));
            plot(xr/1e6, max(f_pem,0)*1e6, '--', 'Color', method_colors{2}, 'LineWidth', 2);

            mu_cq = mean_Q_cum(l); sig_cq = sqrt(max(var_Q_cum(l),0)); gam_cq = skew_Q_cum(l);
            z_cq = (xr - mu_cq) / max(sig_cq, eps);
            f_cum = normpdf(z_cq)/max(sig_cq,eps) .* (1 + gam_cq/6*(z_cq.^3 - 3*z_cq));
            plot(xr/1e6, max(f_cum,0)*1e6, ':', 'Color', method_colors{3}, 'LineWidth', 2);

            xlabel('无功功率 (MVar)', 'FontName',fontName_cn, 'FontSize',fontSZ);
            title(sprintf('%s\x2192%s %s相 Q', info_t.from_bus, info_t.to_bus, phaseChar_arr{ph}), ...
                'FontName',fontName_cn, 'FontSize',fontSZ);
        end

        hold off;
        set(gca, 'fontsize',fontSZ, 'fontname',fontName_en, 'box','off');
        ylabel('PDF', 'FontSize',fontSZ);
        legend({'MCS','本文','PEM','CM'}, 'FontName',fontName_cn, 'FontSize',fontSZ-2, 'Location','best');
        saveas(gcf, fullfile(resultDir, sprintf('Fig04_%s.png', fig4_suffixes{si})));
    end
end
fprintf('[图4] 已保存\n');

% === 表1: 全网统计矩误差与耗时对比表（按方法分行）===
t_per_node_exact_ms = t_vmag_exact / nc * 1000;
t_per_node_delta_ms = t_vmag_delta / nc * 1000;
t_per_node_pem_ms   = t_pem_dss / nc * 1000;
t_per_node_cum_ms   = t_cum_dss / nc * 1000;

fprintf('\n');
fprintf('================================================================================================================================\n');
fprintf('  表1: 全网基础物理量统计矩最大误差 (%%) 与单节点平均耗时\n');
fprintf('================================================================================================================================\n');
fprintf('  %-18s | %8s %8s | %8s %8s | %8s %8s | %10s\n', ...
    '方法', 'Vmag均值', 'Vmag方差', 'P均值', 'P方差', 'Q均值', 'Q方差', '耗时(ms/节点)');
fprintf('  -------------------|-------------------|-------------------|-------------------|-----------\n');
fprintf('  %-18s | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %10.2f\n', '本文-精确Laplace', ...
    max(error_Vmag_mean_exact), max(error_Vmag_var_exact), ...
    max(errMP_aqppf), max(errVP_aqppf), ...
    max(errMQ_aqppf), max(errVQ_aqppf), t_per_node_exact_ms);
fprintf('  %-18s | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %10.4f\n', '本文-Delta近似', ...
    max(error_Vmag_mean_approx), max(error_Vmag_var_approx), ...
    max(errMP_aqppf), max(errVP_aqppf), ...
    max(errMQ_aqppf), max(errVQ_aqppf), t_per_node_delta_ms);
fprintf('  %-18s | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %10.2f\n', 'PEM', ...
    max(errVmag_mean_pem), max(errVmag_var_pem), ...
    max(errMP_pem), max(errVP_pem), ...
    max(errMQ_pem), max(errVQ_pem), t_per_node_pem_ms);
fprintf('  %-18s | %8.4f %8.4f | %8.4f %8.4f | %8.4f %8.4f | %10.2f\n', 'CM', ...
    max(errVmag_mean_cum), max(errVmag_var_cum), ...
    max(errMP_cum), max(errVP_cum), ...
    max(errMQ_cum), max(errVQ_cum), t_per_node_cum_ms);
fprintf('================================================================================================================================\n');

% === Fig.5: 序电压PDF对比 (正序/负序/零序, MCS/本文/PEM/CM) ===
fprintf('\n[图5] 序电压PDF对比...\n');

% PEM/CM 序电压计算
mean_Vpos_pem = zeros(n3ph,1); mean_Vneg_pem = zeros(n3ph,1); mean_V0_pem = zeros(n3ph,1);
var_Vpos_pem  = zeros(n3ph,1); var_Vneg_pem  = zeros(n3ph,1); var_V0_pem  = zeros(n3ph,1);
skew_Vpos_pem = zeros(n3ph,1); skew_Vneg_pem = zeros(n3ph,1); skew_V0_pem = zeros(n3ph,1);
Vpos_cum_all  = zeros(n3ph, N); Vneg_cum_all = zeros(n3ph, N); V0_cum_all = zeros(n3ph, N);

for b = 1:n3ph
    idxA = threePhaseBusIdx(b,1); idxB = threePhaseBusIdx(b,2); idxC = threePhaseBusIdx(b,3);

    % PEM: 加权矩
    Vc_abc = [V_pem_center(idxA); V_pem_center(idxB); V_pem_center(idxC)];
    Vs_c = T_seq * Vc_abc;
    vpos_c = abs(Vs_c(2)); vneg_c = abs(Vs_c(3)); v0_c = abs(Vs_c(1));

    vpos_pts = zeros(2*m_rv,1); vneg_pts = zeros(2*m_rv,1); v0_pts = zeros(2*m_rv,1);
    for kk = 1:2*m_rv
        Vk = [V_pem_all(idxA,kk); V_pem_all(idxB,kk); V_pem_all(idxC,kk)];
        Vsk = T_seq * Vk;
        vpos_pts(kk) = abs(Vsk(2)); vneg_pts(kk) = abs(Vsk(3)); v0_pts(kk) = abs(Vsk(1));
    end
    mean_Vpos_pem(b) = w0_pem*vpos_c + sum(w_pem.*vpos_pts);
    mean_Vneg_pem(b) = w0_pem*vneg_c + sum(w_pem.*vneg_pts);
    mean_V0_pem(b)   = w0_pem*v0_c   + sum(w_pem.*v0_pts);
    EV2_pos = w0_pem*vpos_c^2 + sum(w_pem.*vpos_pts.^2);
    EV2_neg = w0_pem*vneg_c^2 + sum(w_pem.*vneg_pts.^2);
    EV2_0   = w0_pem*v0_c^2   + sum(w_pem.*v0_pts.^2);
    var_Vpos_pem(b) = max(EV2_pos - mean_Vpos_pem(b)^2, 0);
    var_Vneg_pem(b) = max(EV2_neg - mean_Vneg_pem(b)^2, 0);
    var_V0_pem(b)   = max(EV2_0   - mean_V0_pem(b)^2,   0);
    EV3_pos = w0_pem*vpos_c^3 + sum(w_pem.*vpos_pts.^3);
    EV3_neg = w0_pem*vneg_c^3 + sum(w_pem.*vneg_pts.^3);
    EV3_0   = w0_pem*v0_c^3   + sum(w_pem.*v0_pts.^3);
    skew_Vpos_pem(b) = (EV3_pos - 3*mean_Vpos_pem(b)*EV2_pos + 2*mean_Vpos_pem(b)^3) / max(var_Vpos_pem(b)^1.5, eps);
    skew_Vneg_pem(b) = (EV3_neg - 3*mean_Vneg_pem(b)*EV2_neg + 2*mean_Vneg_pem(b)^3) / max(var_Vneg_pem(b)^1.5, eps);
    skew_V0_pem(b)   = (EV3_0   - 3*mean_V0_pem(b)*EV2_0     + 2*mean_V0_pem(b)^3)   / max(var_V0_pem(b)^1.5,   eps);

    % CM: 从线性化电压 V_cum_complex 计算序电压
    V_abc_cum = [V_cum_complex(idxA,:); V_cum_complex(idxB,:); V_cum_complex(idxC,:)];
    V_seq_cum = T_seq * V_abc_cum;
    Vpos_cum_all(b,:) = abs(V_seq_cum(2,:));
    Vneg_cum_all(b,:) = abs(V_seq_cum(3,:));
    V0_cum_all(b,:)   = abs(V_seq_cum(1,:));
end

show_seqv_buses = {'634', '671'};
show_3ph_fig5 = [];
for idx_tmp = 1:length(show_seqv_buses)
    match_tmp = find(cellfun(@(x) strcmpi(x, show_seqv_buses{idx_tmp}), threePhaseBuses));
    assert(~isempty(match_tmp), '未找到三相母线 %s', show_seqv_buses{idx_tmp});
    show_3ph_fig5 = [show_3ph_fig5, match_tmp(1)];
end
n_seq_samp_fig5 = 1e5;
seq_names = {'正序 |V+|', '负序 |V-|', '零序 |V0|'};
seq_mc_data = {Vmag_pos_mc, Vmag_neg_mc, Vmag_zero_mc};
seq_lumda = {Vpos2_lumda, Vneg2_lumda, V0sq_lumda};
seq_delta = {Vpos2_delta, Vneg2_delta, V0sq_delta};
seq_const = {Vpos2_const, Vneg2_const, V0sq_const};
pem_seq_mean = {mean_Vpos_pem, mean_Vneg_pem, mean_V0_pem};
pem_seq_var  = {var_Vpos_pem,  var_Vneg_pem,  var_V0_pem};
pem_seq_skew = {skew_Vpos_pem, skew_Vneg_pem, skew_V0_pem};
cum_seq_all  = {Vpos_cum_all,  Vneg_cum_all,  V0_cum_all};
cum_seq_mean = cellfun(@(x) mean(x,2), cum_seq_all, 'UniformOutput', false);
cum_seq_var  = cellfun(@(x) var(x,0,2), cum_seq_all, 'UniformOutput', false);
cum_seq_skew = cell(1,3);
for si_tmp = 1:3
    s_tmp = zeros(n3ph,1);
    for b_tmp = 1:n3ph, s_tmp(b_tmp) = skewness(cum_seq_all{si_tmp}(b_tmp,:)'); end
    cum_seq_skew{si_tmp} = s_tmp;
end

seq_fig_suffix = {'pos', 'neg', 'zero'};
for bi = 1:length(show_3ph_fig5)
    b = show_3ph_fig5(bi);
    Vb_seq = V_base_node(threePhaseBusIdx(b,1));
    busLabel_fig5 = busLabels_of_node{threePhaseBusIdx(b,1)};
    for si = 1:3
        figure('Units','pixels','Position',[100+80*si, 100+80*si, 700, 500]);
        mc_data = seq_mc_data{si}(:, b);
        mc_data_pu = mc_data / Vb_seq;
        histogram(mc_data_pu, 60, 'Normalization','pdf','FaceAlpha',0.3, ...
            'FaceColor',[0.07 0.44 0.75], 'EdgeColor','none');
        hold on;

        % 本文: GNCS采样 -> ksdensity
        smp_plot = [];
        for k = 1:gmm_num
            nk = round(n_seq_samp_fig5*com_pi(k)); if nk==0, continue; end
            lk = seq_lumda{si}{b,k}; dk = seq_delta{si}{b,k};
            Z = randn(length(lk), nk);
            sq = seq_const{si}{b,k} + sum(lk.*(Z+sqrt(dk)).^2,1);
            smp_plot = [smp_plot, sqrt(max(sq,0))];
        end
        smp_plot_pu = smp_plot / Vb_seq;
        [pdf_aqppf, xr_aqppf] = ksdensity(smp_plot_pu, 'NumPoints', 300);
        plot(xr_aqppf, pdf_aqppf, '-', 'Color', method_colors{1}, 'LineWidth', 4);

        % PEM: Gram-Charlier
        mu_seq = pem_seq_mean{si}(b)/Vb_seq;
        sig_seq = sqrt(max(pem_seq_var{si}(b),0))/Vb_seq;
        xr_seq = linspace(min(mc_data_pu)*0.95, max(mc_data_pu)*1.05, 300);
        z_seq = (xr_seq - mu_seq) / max(sig_seq, eps);
        H3_seq = z_seq.^3 - 3*z_seq;
        f_pem_seq = normpdf(z_seq)/max(sig_seq,eps) .* (1 + pem_seq_skew{si}(b)/6*H3_seq);
        f_pem_seq = max(f_pem_seq, 0);
        plot(xr_seq, f_pem_seq, '--', 'Color', method_colors{2}, 'LineWidth', 4);

        % CM: Gram-Charlier
        mu_cum_seq = cum_seq_mean{si}(b)/Vb_seq;
        sig_cum_seq = sqrt(max(cum_seq_var{si}(b),0))/Vb_seq;
        z_cum_seq = (xr_seq - mu_cum_seq) / max(sig_cum_seq, eps);
        H3_cum_seq = z_cum_seq.^3 - 3*z_cum_seq;
        f_cum_seq = normpdf(z_cum_seq)/max(sig_cum_seq,eps) .* (1 + cum_seq_skew{si}(b)/6*H3_cum_seq);
        f_cum_seq = max(f_cum_seq, 0);
        plot(xr_seq, f_cum_seq, ':', 'Color', method_colors{3}, 'LineWidth', 4);

        hold off;
        set(gca,'fontsize',fontFig,'fontname',fontName_en,'box','off');
        if si == 1 && contains(busLabel_fig5, '634', 'IgnoreCase', true)
            xlim([0.999, 1.015]);
        elseif si == 1 && contains(busLabel_fig5, '671', 'IgnoreCase', true)
            xlim([1.01, 1.025]);
        end
        xlabel('电压幅值 (p.u.)','FontName',fontName_cn,'FontSize',fontFig);
        ylabel('PDF','FontSize',fontFig);
        title(sprintf('%s %s', busLabel_fig5, seq_names{si}),'FontName',fontName_cn,'FontSize',fontFig);
        legend({'MCS','本文','PEM','CM'},'FontName',fontName_cn,'FontSize',fontFig-2,'Location','best');
        saveas(gcf, fullfile(resultDir, sprintf('Fig05_%s_SeqV_%s.png', busLabel_fig5, seq_fig_suffix{si})));
    end
end
fprintf('[图5] 已保存 (3张独立图/母线)\n');

% === Fig.6: VUF PDF对比 (MCS/本文/PEM/CM)，所有三相节点各一图 ===
fprintf('\n[图6] VUF PDF对比，共 %d 个三相节点...\n', n3ph);
fontVUF = 22;
n_pdf_vuf = 300;

for idx = 1:n3ph
    b = idx;
    busLabel_b = busLabels_of_node{threePhaseBusIdx(b,1)};
    if any(strcmpi(busLabel_b, {'PCC','Reg','650'})), continue; end

    y_lo = max(0.001, min(VUF_mc(:,b)) * 0.8);
    y_hi = max(VUF_mc(:,b)) * 1.2;
    y_pts = linspace(y_lo, y_hi, n_pdf_vuf);

    figure('Units','pixels','Position',[100+60*idx, 100+60*idx, 700, 500]);
    histogram(VUF_mc(:,b), 80, 'Normalization','pdf', 'FaceAlpha',0.3, ...
        'FaceColor',[0.07 0.44 0.75], 'EdgeColor','none');
    hold on;

    % 本文: 解析PDF
    f_vuf_exact = zeros(1, n_pdf_vuf);
    for k = 1:gmm_num
        f_k = vuf_pdf_eq48(y_pts(:), VUF_lam_pos(b,k), VUF_lam_neg(b,k), ...
            VUF_rho(b,k), VUF_sig_ratio(b,k));
        f_vuf_exact = f_vuf_exact + com_pi(k) * f_k(:)';
    end
    plot(y_pts, f_vuf_exact, '-', 'Color', method_colors{1}, 'LineWidth', 4);

    % PEM: Gram-Charlier
    mu_p = mean_VUF_pem(b); sig_p = sqrt(max(var_VUF_pem(b), eps)); gam_p = skew_VUF_pem(b);
    z_p = (y_pts - mu_p) / sig_p;
    f_pem = normpdf(z_p)/sig_p .* (1 + gam_p/6*(z_p.^3 - 3*z_p));
    f_pem = max(f_pem, 0);
    plot(y_pts, f_pem, '--', 'Color', method_colors{2}, 'LineWidth', 4);

    % CM: Gram-Charlier
    mu_cv = mean_VUF_cum(b); sig_cv = sqrt(max(var_VUF_cum(b), eps)); gam_cv = skew_VUF_cum(b);
    z_cv = (y_pts - mu_cv) / sig_cv;
    f_cum = normpdf(z_cv)/sig_cv .* (1 + gam_cv/6*(z_cv.^3 - 3*z_cv));
    f_cum = max(f_cum, 0);
    plot(y_pts, f_cum, ':', 'Color', method_colors{3}, 'LineWidth', 4);

    hold off;
    set(gca, 'FontSize', fontVUF, 'FontName', fontName_en, 'Box', 'off');
    xlim([0 max(y_hi, 1)]);
    xlabel('VUF (%)', 'FontName', fontName_cn, 'FontSize', fontVUF);
    ylabel('PDF', 'FontSize', fontVUF, 'FontName', fontName_en);
    legend({'MCS','本文','PEM','CM'}, 'FontName', fontName_cn, 'FontSize', fontVUF, 'Location', 'best');
    title(sprintf('%s VUF概率分布对比', busLabel_b), 'FontName', fontName_cn, 'FontSize', fontVUF);
    saveas(gcf, fullfile(resultDir, sprintf('Fig06_VUF_PDF_%s.png', busLabel_b)));
end
fprintf('[图6] 已保存\n');

% === Fig.6b: VUF CDF对比 (式51 解析 / MCS ecdf / PEM / CM)，所有三相节点各一图 ===
fprintf('\n[图6b] VUF CDF对比 (式51 解析)，共 %d 个三相节点...\n', n3ph);
n_cdf_vuf = 200;

for idx = 1:n3ph
    b = idx;
    busLabel_b = busLabels_of_node{threePhaseBusIdx(b,1)};
    if any(strcmpi(busLabel_b, {'PCC','Reg','650'})), continue; end

    ups_max = max(VUF_mc(:,b)) * 1.2 / 100;
    ups_pts = linspace(0, ups_max, n_cdf_vuf);

    % 本文: 式(51) 参数扫描 — F_VUF(υ) = F_{Z(υ)}(0), Z(υ)=V'(M⁻−υ²M⁺)V
    F_vuf_analytical = zeros(1, n_cdf_vuf);
    fprintf('  [式51] 扫描 %d 个υ值 (母线 %s)...\n', n_cdf_vuf, busLabel_b);
    for ui = 1:n_cdf_vuf
        ups = ups_pts(ui);
        if ups < 1e-12
            F_vuf_analytical(ui) = 0;
            continue;
        end
        M_diff = M_Vneg2{b} - ups^2 * M_Vpos2{b};
        F_k_vals = zeros(gmm_num, 1);
        for k = 1:gmm_num
            mu_k = gmm_fit.mu(k,:)';
            BMB = B_GMM{k}' * M_diff * B_GMM{k};
            BMB_sym = 0.5*(BMB + BMB');
            [Q_d, D_d] = eig(BMB_sym);
            lv = diag(D_d);
            a_d = 2 * B_GMM{k}' * M_diff * mu_k;
            bv = Q_d' * a_d;
            alpha_d = mu_k' * M_diff * mu_k;
            tol_d = 1e-10 * max(abs(lv));
            inz = abs(lv) > max(tol_d, 1e-30);
            if ~any(inz)
                F_k_vals(k) = double(alpha_d <= 0);
                continue;
            end
            w_row = lv(inz)';
            k_row = ones(1, sum(inz));
            cv = zeros(size(lv));
            cv(inz) = bv(inz) ./ (2*lv(inz));
            dv = cv.^2;
            lambda_row = dv(inz)';
            m_const = alpha_d - sum(lv(inz).*dv(inz));
            try
                F_k_vals(k) = gx2cdf(0, w_row, k_row, lambda_row, m_const, 0);
            catch
                F_k_vals(k) = gx2cdf_davies(0, w_row, k_row, lambda_row, m_const, 0);
            end
            F_k_vals(k) = max(0, min(1, F_k_vals(k)));
        end
        F_vuf_analytical(ui) = com_pi * F_k_vals;
    end

    figure('Units','pixels','Position',[200+60*idx, 200+60*idx, 700, 500]);
    hold on;

    % 本文: 式(51) 解析CDF
    plot(ups_pts*100, F_vuf_analytical, '-*', 'Color', method_colors{1}, 'LineWidth', 4, ...
        'MarkerIndices', 1:10:length(F_vuf_analytical));

    % PEM: Gram-Charlier CDF = Phi(z) - phi(z)*[gamma/6*(z^2-1)]
    mu_p = mean_VUF_pem(b); sig_p = sqrt(max(var_VUF_pem(b), eps)); gam_p = skew_VUF_pem(b);
    x_cdf = linspace(0, max(VUF_mc(:,b))*1.2, 300);
    z_p = (x_cdf - mu_p) / sig_p;
    F_pem = normcdf(z_p) - normpdf(z_p) .* (gam_p/6*(z_p.^2 - 1));
    F_pem = max(0, min(1, F_pem));
    plot(x_cdf, F_pem, '--', 'Color', method_colors{2}, 'LineWidth', 4);

    % CM: ecdf from linearized samples (针对当前母线b重新计算序电压)
    idxA_b = threePhaseBusIdx(b,1); idxB_b = threePhaseBusIdx(b,2); idxC_b = threePhaseBusIdx(b,3);
    V_abc_cum_b = [V_cum_complex(idxA_b,:); V_cum_complex(idxB_b,:); V_cum_complex(idxC_b,:)];
    V_seq_cum_b = T_seq * V_abc_cum_b;
    VUF_cum_b = abs(V_seq_cum_b(3,:)) ./ abs(V_seq_cum_b(2,:)) * 100;
    [F_cum_cdf, x_cum_cdf] = ecdf(VUF_cum_b);
    plot(x_cum_cdf, F_cum_cdf, ':', 'Color', method_colors{3}, 'LineWidth', 4);

    % MCS ecdf — 最后画，置于最上层
    [F_mc, x_mc] = ecdf(VUF_mc(:,b));
    stairs(x_mc, F_mc, '-', 'Color', [0 0 0], 'LineWidth', 2);

    hold off;
    set(gca, 'FontSize', fontVUF, 'FontName', fontName_en, 'Box', 'off');
    xlim([0 max(ups_max*100*1.1, 1)]);
    xlabel('VUF (%)', 'FontName', fontName_cn, 'FontSize', fontVUF);
    ylabel('CDF', 'FontSize', fontVUF, 'FontName', fontName_en);
    legend({'本文','PEM','CM','MCS'}, 'FontName', fontName_cn, 'FontSize', fontVUF-2, 'Location', 'best');
    title(sprintf('%s VUF累积分布对比', busLabel_b), 'FontName', fontName_cn, 'FontSize', fontVUF);
    saveas(gcf, fullfile(resultDir, sprintf('Fig06b_VUF_CDF_%s.png', busLabel_b)));
end
fprintf('[图6b] 已保存\n');

% === 表2: VUF统计矩详细对比表 ===
fprintf('\n');
fprintf('============================================================================================================\n');
fprintf('  表2: VUF统计矩详细对比\n');
fprintf('============================================================================================================\n');
fprintf('  %-8s | %10s %10s %10s %10s %10s | %10s %10s %10s %10s %10s\n', ...
    '母线', 'MCS均值', '精确积分', 'Delta', 'PEM', 'CM', ...
    'MCS方差', '精确积分', 'Delta', 'PEM', 'CM');
fprintf('  ---------|-----------------------------------------------------|-----------------------------------------------------\n');
[~, vuf_sort_tab2] = sort(threePhaseBusIdx(:,1));
for ii = 1:n3ph
    b = vuf_sort_tab2(ii);
    fprintf('  %-8s | %10.4f %10.4f %10.4f %10.4f %10.4f | %10.6f %10.6f %10.6f %10.6f %10.6f\n', ...
        busLabels_3ph_short{ii}, ...
        mean_VUF_mc(b), mean_VUF_exact(b), mean_VUF_approx(b), mean_VUF_pem(b), mean_VUF_cum(b), ...
        var_VUF_mc(b), var_VUF_exact(b), var_VUF_approx(b), var_VUF_pem(b), var_VUF_cum(b));
end
fprintf('============================================================================================================\n');

%% ========================================================================
%  汇总输出
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  AQPPF_3ph 计算完成\n');
fprintf('%s\n', repmat('=',1,70));
fprintf('  算例: %s\n', caseName);
fprintf('  计算节点数: %d\n', nc);
fprintf('  三相节点数: %d\n', n3ph);
fprintf('  GMM分量数: %d\n', gmm_num);
fprintf('  线路功率均值最大误差: P=%.4f%%, Q=%.4f%%\n', ...
    max(error_Line_P_mean), max(error_Line_Q_mean));
fprintf('  线路功率方差最大误差: P=%.4f%%, Q=%.4f%%\n', ...
    max(error_Line_P_var), max(error_Line_Q_var));
fprintf('  电压幅值均值误差(精确/近似): %.4f%% / %.4f%%\n', ...
    max(error_Vmag_mean_exact), max(error_Vmag_mean_approx));
fprintf('  电压幅值方差误差(精确/近似): %.4f%% / %.4f%%\n', ...
    max(error_Vmag_var_exact), max(error_Vmag_var_approx));
fprintf('  序电压幅值均值最大误差: V+=%.4f%%, V-=%.4f%%, V0=%.4f%%\n', ...
    max(err_Vpos_mean), max(err_Vneg_mean), max(err_V0_mean));
fprintf('  VUF均值误差(精确/近似): %.4f%% / %.4f%%\n', ...
    max(error_VUF_mean_exact), max(error_VUF_mean_approx));
fprintf('  VUF方差误差(精确/近似): %.4f%% / %.4f%%\n', ...
    max(error_VUF_var_exact), max(error_VUF_var_approx));
fprintf('\n  --- 对比方法汇总 ---\n');
fprintf('  P均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errMP_aqppf), max(errMP_pem), max(errMP_cum));
fprintf('  P方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVP_aqppf), max(errVP_pem), max(errVP_cum));
fprintf('  VUF均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVUF_mean_aqppf), max(errVUF_mean_pem), max(errVUF_mean_cum));
fprintf('  VUF方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVUF_var_aqppf), max(errVUF_var_pem), max(errVUF_var_cum));
fprintf('  电压幅值均值最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVmag_mean_aqppf), max(errVmag_mean_pem), max(errVmag_mean_cum));
fprintf('  电压幅值方差最大误差: 本文=%.2f%%, PEM=%.2f%%, CM=%.2f%%\n', ...
    max(errVmag_var_aqppf), max(errVmag_var_pem), max(errVmag_var_cum));
fprintf('  计算时间: 本文=%.1fs, PEM=%.1fs, CM=%.1fs\n', ...
    t3, t_pem_dss, t_cum_dss);
fprintf('  结果保存目录: %s\n', resultDir);
fprintf('%s\n', repmat('=',1,70));

% 保存工作空间
wsFile = fullfile(resultDir, 'workspace.mat');
if isfile(wsFile), delete(wsFile); end
save(wsFile, '-v7.3');
fprintf('[保存] 工作空间已保存\n');
