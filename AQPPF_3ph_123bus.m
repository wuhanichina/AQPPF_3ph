% AQPPF_3PH_123BUS.M
% IEEE 123节点系统的解析二次型概率潮流计算
% 含PEM对比方法
% 产物: Fig.7/7b (全网误差散点), Fig.8 (K值敏感性), Tab.3 (耗时拆解)
%
% 版本: v1.0
% 日期: 2026-02-25

dbstop if error; clc; clear; close all;

%% ========================================================================
%  段落1: 初始化
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  AQPPF_3ph_123bus - IEEE 123节点系统\n');
fprintf('%s\n\n', repmat('=',1,70));

projectDir = fileparts(mfilename('fullpath'));
addpath(fullfile(projectDir, 'gncs'));
addpath(projectDir);

load(fullfile(projectDir, 'data', 'Load_Data.mat'));

caseName = 'IEEE123';
gmm_num  = 20;
sample_num = 10000;
dssFile = 'C:\Program Files\OpenDSS\IEEETestCases\123Bus\IEEE123Master.dss';
assert(isfile(dssFile), '[依赖缺失] OpenDSS算例文件不存在: %s', dssFile);

resultDir = fullfile(projectDir, 'result', caseName);
cacheDir  = fullfile(projectDir, 'cache');
if ~isfolder(resultDir), mkdir(resultDir); end
if ~isfolder(cacheDir),  mkdir(cacheDir);  end

fprintf('[配置] 算例: %s\n', caseName);
fprintf('[配置] GMM分量数 K=%d, 采样数=%d\n', gmm_num, sample_num);

%% ========================================================================
%  段落2: 蒙特卡洛模拟（OpenDSS三相潮流，带缓存）
%  ========================================================================
N = height(load_data);
cacheFile = fullfile(cacheDir, sprintf('MC_3ph_%s_T%d.mat', caseName, N));

if isfile(cacheFile)
    fprintf('\n[缓存] 发现缓存文件: %s\n', cacheFile);
    cacheData = load(cacheFile);
    V_mc_all     = cacheData.V_mc_all;
    nodeNames    = cacheData.nodeNames;
    nc           = cacheData.nc;
    Y3ph         = cacheData.Y3ph;
    nodeOrderY   = cacheData.nodeOrderY;
    lineNames    = cacheData.lineNames;
    nLine        = cacheData.nLine;
    busNames     = cacheData.busNames;
    t1 = 0;
    fprintf('[缓存] 加载完成: %d 节点, %d 样本\n', nc, N);

    % 重建负荷信息（PEM对比需要）
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
    fprintf('\n[MC] 运行蒙特卡洛模拟 (N=%d)...\n', N);
    DSSObj = actxserver('OpenDSSEngine.DSS');
    assert(DSSObj.Start(0), '[OpenDSS] 启动失败');
    DSSText = DSSObj.Text;

    DSSText.Command = ['Compile "', dssFile, '"'];
    DSSCircuit = DSSObj.ActiveCircuit;

    iLoad = DSSCircuit.Loads;
    nLoad = iLoad.Count;
    loadNames = cell(nLoad, 1);
    loadKW_base  = zeros(nLoad, 1);
    loadKvar_base = zeros(nLoad, 1);
    idx = iLoad.First; li = 1;
    while idx > 0
        loadNames{li}    = iLoad.Name;
        loadKW_base(li)  = iLoad.kW;
        loadKvar_base(li) = iLoad.kvar;
        li = li + 1; idx = iLoad.Next;
    end

    nodeNames = DSSCircuit.AllNodeNames;
    nc = length(nodeNames);
    busNames  = DSSCircuit.AllBusNames;

    iLine = DSSCircuit.Lines;
    nLine = iLine.Count;
    lineNames = cell(nLine, 1);
    idx = iLine.First; li = 1;
    while idx > 0
        lineNames{li} = iLine.Name;
        li = li + 1; idx = iLine.Next;
    end

    sysY = DSSCircuit.SystemY;
    Y_vec = sysY(1:2:end) + 1i * sysY(2:2:end);
    nY = round(sqrt(length(Y_vec)));
    Y3ph = reshape(Y_vec, [nY, nY]);
    nodeOrderY = DSSCircuit.YNodeOrder;

    fprintf('[OpenDSS] nc=%d, nLine=%d\n', nc, nLine);

    load_sequence = mod(0:nLoad-1, 24) + 1;
    V_mc_all = zeros(N, nc);

    tic;
    h = waitbar(0, sprintf('MC: 0/%d', N), 'Name', 'IEEE123 MC');
    for t = 1:N
        DSSText.Command = ['Compile "', dssFile, '"'];
        DSSCircuit = DSSObj.ActiveCircuit;

        iLoad = DSSCircuit.Loads;
        idx = iLoad.First; li = 1;
        while idx > 0 && li <= nLoad
            scaleFactor = load_data{t, load_sequence(li)};
            iLoad.kW   = loadKW_base(li) * scaleFactor;
            iLoad.kvar = loadKvar_base(li) * scaleFactor;
            li = li + 1; idx = iLoad.Next;
        end

        DSSSolution = DSSCircuit.Solution;
        DSSSolution.Solve;

        allVolts = DSSCircuit.AllBusVolts;
        V_mc_all(t, :) = allVolts(1:2:end) + 1i * allVolts(2:2:end);

        if mod(t, 200) == 0
            waitbar(t/N, h, sprintf('MC: %d/%d', t, N));
        end
    end
    close(h);
    t1 = toc;

    save(cacheFile, 'V_mc_all', 'nodeNames', 'nc', 'Y3ph', 'nodeOrderY', ...
         'lineNames', 'nLine', 'busNames', '-v7.3');
    DSSObj.delete; clear DSSObj;
end

if t1 > 0
    fprintf('[MC] MCS耗时: %.2f s\n', t1);
    t_mc = t1;
else
    fprintf('[MC] 已从缓存加载\n');
    t_mc = 0;
end
fprintf('[MC] N=%d, nc=%d\n', N, nc);

%% ========================================================================
%  段落2b: 全局绘图设置与节点编号
%  ========================================================================
fontSZ = 14; fontSZ_title = 16;
fontName_cn = '宋体'; fontName_en = 'Times New Roman';
phaseChar_arr = {'A','B','C'};

nodeLabels = cell(1, nc);
busLabels_of_node = cell(1, nc);
phaseLabels_of_node = zeros(1, nc);
for i = 1:nc
    parts = split(nodeNames{i}, '.');
    busRaw = lower(parts{1});
    phNum  = str2double(parts{2});
    busLabel = upper(busRaw);
    busLabels_of_node{i}    = busLabel;
    phaseLabels_of_node(i)  = phNum;
    nodeLabels{i} = sprintf('%s-%s', busLabel, phaseChar_arr{phNum});
end

%% ========================================================================
%  段落3: 电压实同构
%  ========================================================================
v_vec = [real(V_mc_all), imag(V_mc_all)];  % [N x 2*nc]
V_base_node = mean(abs(V_mc_all), 1);  % [1 x nc]
fprintf('[实同构] 完成: [%d x %d]\n', size(v_vec));

%% ========================================================================
%  段落4: GMM拟合
%  ========================================================================
fprintf('\n[GMM] 拟合 %d 分量GMM...\n', gmm_num);
tic;
gmm_fit = fitgmdist(v_vec, gmm_num, 'RegularizationValue', 1e-6, ...
    'Options', statset('MaxIter', 500, 'TolFun', 1e-8), ...
    'CovarianceType', 'full', 'SharedCovariance', false);
t_gmm = toc;
com_pi = gmm_fit.ComponentProportion;
fprintf('[GMM] 完成, 耗时: %.2f s\n', t_gmm);

%% ========================================================================
%  段落5: 二次型参数矩阵
%  ========================================================================
fprintf('\n[二次型] 构建矩阵...\n');

B_GMM = cell(gmm_num, 1);
for k = 1:gmm_num
    Sigma_k = gmm_fit.Sigma(:,:,k);
    [V_eig, D_eig] = eig(Sigma_k);
    d = diag(D_eig); d(d < 0) = 0;
    B_GMM{k} = V_eig * diag(sqrt(d));
end

nodeOrderY_lower = lower(nodeOrderY);
nodeNames_lower  = lower(nodeNames);
permY2V = zeros(nc, 1);
for i = 1:nc
    idx_match = find(strcmp(nodeOrderY_lower{i}, nodeNames_lower));
    assert(~isempty(idx_match), '节点 %s 未找到', nodeOrderY{i});
    permY2V(i) = idx_match;
end
permV2Y = zeros(nc, 1);
permV2Y(permY2V) = 1:nc;
Y3ph_ordered = Y3ph(permV2Y, permV2Y);

[Node_P_matrix, Node_Q_matrix] = get_3ph_node_power_matrix(Y3ph_ordered, nc);

branchPairs_all = [];
for i = 1:nc
    for j = i+1:nc
        if abs(Y3ph_ordered(i,j)) > 1e-10
            branchPairs_all = [branchPairs_all; i, j];
        end
    end
end

Y_adm_threshold = 50;
isExcluded = false(size(branchPairs_all, 1), 1);
for l = 1:size(branchPairs_all, 1)
    ni = branchPairs_all(l,1); nj = branchPairs_all(l,2);
    Vratio = V_base_node(ni) / V_base_node(nj);
    Ymag = abs(Y3ph_ordered(ni, nj));
    if Vratio > 1.5 || Vratio < 1/1.5
        isExcluded(l) = true;
    elseif Ymag > Y_adm_threshold
        isExcluded(l) = true;
    end
end
branchPairs = branchPairs_all(~isExcluded, :);
nBranch = size(branchPairs, 1);
fprintf('[线路识别] 保留线路: %d\n', nBranch);

[Line_P_matrix, Line_Q_matrix] = get_3ph_line_power_matrix(Y3ph_ordered, branchPairs, nc);
fprintf('[二次型] 线路功率矩阵: %d 个\n', nBranch);

%% ========================================================================
%  段落6: GNCS参数计算
%  ========================================================================
fprintf('\n[GNCS] 计算线路功率参数...\n');

Line_P_lumda = cell(nBranch, gmm_num);
Line_P_delta = cell(nBranch, gmm_num);
Line_P_const = cell(nBranch, gmm_num);
Line_Q_lumda = cell(nBranch, gmm_num);
Line_Q_delta = cell(nBranch, gmm_num);
Line_Q_const = cell(nBranch, gmm_num);

tic;
h = waitbar(0, sprintf('GNCS: 0/%d', nBranch), 'Name', 'GNCS参数');
for l = 1:nBranch
    for k = 1:gmm_num
        mu_k = gmm_fit.mu(k,:)';
        for pq = 1:2
            if pq == 1, M = Line_P_matrix{l}; else, M = Line_Q_matrix{l}; end
            BMB = B_GMM{k}' * M * B_GMM{k};
            BMB_sym = 0.5*(BMB + BMB');
            [Q_m, D_m] = eig(BMB_sym);
            lm = diag(D_m);
            am = 2 * B_GMM{k}' * M * mu_k;
            bm = Q_m' * am;
            alpha_m = mu_k' * M * mu_k;
            tol = 1e-12 * max(abs(lm));
            inz = abs(lm) > tol;
            cm = zeros(size(lm)); cm(inz) = bm(inz)./(2*lm(inz));
            dm = cm.^2;
            if pq == 1
                Line_P_lumda{l,k} = lm; Line_P_delta{l,k} = dm;
                Line_P_const{l,k} = alpha_m - sum(lm(inz).*dm(inz));
            else
                Line_Q_lumda{l,k} = lm; Line_Q_delta{l,k} = dm;
                Line_Q_const{l,k} = alpha_m - sum(lm(inz).*dm(inz));
            end
        end
    end
    if mod(l, 20) == 0
        waitbar(l/nBranch, h, sprintf('GNCS: %d/%d', l, nBranch));
    end
end
close(h);
t_gncs = toc;
fprintf('[GNCS] 耗时: %.2f s\n', t_gncs);

%% ========================================================================
%  段落7: 线路功率统计矩与误差
%  ========================================================================
fprintf('\n[统计矩] 计算线路功率...\n');

samePhaseIdx = [];
for l = 1:nBranch
    parts_i = split(nodeNames{branchPairs(l,1)}, '.');
    parts_j = split(nodeNames{branchPairs(l,2)}, '.');
    if strcmp(parts_i{2}, parts_j{2})
        samePhaseIdx = [samePhaseIdx; l];
    end
end
nSP = length(samePhaseIdx);

V_mc_T = v_vec';
Line_Power_P_MC = cell(nBranch, 1);
Line_Power_Q_MC = cell(nBranch, 1);
for l = 1:nBranch
    Line_Power_P_MC{l} = dot(V_mc_T, Line_P_matrix{l} * V_mc_T)';
    Line_Power_Q_MC{l} = dot(V_mc_T, Line_Q_matrix{l} * V_mc_T)';
end

tic;
mean_Line_P_com = zeros(nBranch, gmm_num);
mean_Line_Q_com = zeros(nBranch, gmm_num);
var_Line_P_com  = zeros(nBranch, gmm_num);
var_Line_Q_com  = zeros(nBranch, gmm_num);

for l = 1:nBranch
    for k = 1:gmm_num
        mean_Line_P_com(l,k) = trace(Line_P_matrix{l} * gmm_fit.Sigma(:,:,k)) ...
            + gmm_fit.mu(k,:) * Line_P_matrix{l} * gmm_fit.mu(k,:)';
        mean_Line_Q_com(l,k) = trace(Line_Q_matrix{l} * gmm_fit.Sigma(:,:,k)) ...
            + gmm_fit.mu(k,:) * Line_Q_matrix{l} * gmm_fit.mu(k,:)';
        var_Line_P_com(l,k) = 2*trace((Line_P_matrix{l}*gmm_fit.Sigma(:,:,k))^2) ...
            + 4*gmm_fit.mu(k,:)*Line_P_matrix{l}*gmm_fit.Sigma(:,:,k)*Line_P_matrix{l}*gmm_fit.mu(k,:)';
        var_Line_Q_com(l,k) = 2*trace((Line_Q_matrix{l}*gmm_fit.Sigma(:,:,k))^2) ...
            + 4*gmm_fit.mu(k,:)*Line_Q_matrix{l}*gmm_fit.Sigma(:,:,k)*Line_Q_matrix{l}*gmm_fit.mu(k,:)';
    end
end

mean_Line_P_this = mean_Line_P_com * com_pi';
mean_Line_Q_this = mean_Line_Q_com * com_pi';
var_Line_P_this  = var_Line_P_com * com_pi' + ...
    sum(com_pi .* (mean_Line_P_com - mean_Line_P_this).^2, 2);
var_Line_Q_this  = var_Line_Q_com * com_pi' + ...
    sum(com_pi .* (mean_Line_Q_com - mean_Line_Q_this).^2, 2);

mean_Line_P_mc = cellfun(@mean, Line_Power_P_MC);
mean_Line_Q_mc = cellfun(@mean, Line_Power_Q_MC);
var_Line_P_mc  = cellfun(@var,  Line_Power_P_MC);
var_Line_Q_mc  = cellfun(@var,  Line_Power_Q_MC);

mean_P_sp = mean_Line_P_mc(samePhaseIdx);
thr_P = 0.01 * median(abs(mean_P_sp));
error_Line_P_mean = abs(mean_P_sp - mean_Line_P_this(samePhaseIdx)) ...
    ./ max(abs(mean_P_sp), thr_P) * 100;
error_Line_P_mean(~isfinite(error_Line_P_mean)) = 0;
t_stat = toc;

fprintf('[统计矩] P均值最大误差: %.4f%%\n', max(error_Line_P_mean));

%% ========================================================================
%  段落9: 三相电压幅值分析
%  ========================================================================
fprintf('\n[电压幅值] 计算...\n');

M_Vmag2 = cell(nc, 1);
for i = 1:nc
    M_Vmag2{i} = zeros(2*nc, 2*nc);
    M_Vmag2{i}(i, i) = 1;
    M_Vmag2{i}(nc+i, nc+i) = 1;
end

Vmag2_lumda = cell(nc, gmm_num);
Vmag2_delta = cell(nc, gmm_num);
Vmag2_const = cell(nc, gmm_num);

tic;
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
        cv = zeros(size(lv)); cv(inz) = bv(inz)./(2*lv(inz));
        dv = cv.^2;
        Vmag2_lumda{i,k} = lv(inz);
        Vmag2_delta{i,k} = dv(inz);
        Vmag2_const{i,k} = alpha_v - sum(lv(inz).*dv(inz));
    end
end
t_vmag_gncs = toc;

Vmag_mc = zeros(N, nc);
for i = 1:nc
    Vmag_mc(:,i) = sqrt(v_vec(:,i).^2 + v_vec(:,nc+i).^2);
end
mean_Vmag_mc = mean(Vmag_mc, 1)';
var_Vmag_mc  = var(Vmag_mc, 0, 1)';

% 精确法: GNCS采样 -> sqrt -> 统计矩
n_gncs_samples = 1e5;
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
        c_vals = sqrt(Vmag2_delta{i,k});
        Vsq_samples = Vmag2_const{i,k} + sum(Vmag2_lumda{i,k} .* (Z + c_vals).^2, 1);
        Vmag_gncs_samples = [Vmag_gncs_samples, sqrt(max(Vsq_samples, 0))];
    end
    mean_Vmag_exact(i) = mean(Vmag_gncs_samples);
    var_Vmag_exact(i)  = var(Vmag_gncs_samples);
end
t_vmag_exact = toc(tic_vmag_exact);
fprintf('[电压幅值] 精确法(GNCS采样)完成, 耗时: %.2f s\n', t_vmag_exact);

% Delta近似法: E[|V|] ≈ sqrt(E[|V|²]),  Var(|V|) ≈ Var(|V|²)/(4·E[|V|²])
% See §3.3, Eq.(42)-(46b) in derivation.md
tic_vmag_delta = tic;
mean_Vmag2_ana = zeros(nc, 1);
var_Vmag2_ana  = zeros(nc, 1);
for i = 1:nc
    mean_W_k = zeros(gmm_num, 1);
    var_W_k  = zeros(gmm_num, 1);
    for k = 1:gmm_num
        mean_W_k(k) = Vmag2_const{i,k} + sum(Vmag2_lumda{i,k} .* (1 + Vmag2_delta{i,k}));
        var_W_k(k)  = 2 * sum(Vmag2_lumda{i,k}.^2 .* (1 + 2*Vmag2_delta{i,k}));
    end
    mean_Vmag2_ana(i) = mean_W_k' * com_pi';
    var_Vmag2_ana(i)  = var_W_k' * com_pi' + mean_W_k.^2' * com_pi' - mean_Vmag2_ana(i)^2;
end
mean_Vmag_delta = sqrt(mean_Vmag2_ana);
var_Vmag_delta  = var_Vmag2_ana ./ (4 * mean_Vmag2_ana);
t_vmag_delta = toc(tic_vmag_delta);
fprintf('[电压幅值] Delta近似法完成, 耗时: %.4f s\n', t_vmag_delta);

thr_vmag_mean = 0.01 * median(abs(mean_Vmag_mc));
thr_vmag_var  = 0.01 * median(abs(var_Vmag_mc));
error_Vmag_mean = abs((mean_Vmag_mc - mean_Vmag_exact) ./ max(abs(mean_Vmag_mc), thr_vmag_mean)) * 100;
error_Vmag_var  = abs((var_Vmag_mc - var_Vmag_exact)   ./ max(abs(var_Vmag_mc),  thr_vmag_var))  * 100;
error_Vmag_mean_delta = abs((mean_Vmag_mc - mean_Vmag_delta) ./ max(abs(mean_Vmag_mc), thr_vmag_mean)) * 100;
error_Vmag_var_delta  = abs((var_Vmag_mc - var_Vmag_delta)   ./ max(abs(var_Vmag_mc),  thr_vmag_var))  * 100;

fprintf('[电压幅值] 精确法: 均值最大误差 %.4f%%, 方差最大误差 %.4f%%\n', ...
    max(error_Vmag_mean), max(error_Vmag_var));
fprintf('[电压幅值] Delta法: 均值最大误差 %.4f%%, 方差最大误差 %.4f%%\n', ...
    max(error_Vmag_mean_delta), max(error_Vmag_var_delta));

%% ========================================================================
%  段落9b: 节点注入功率统计矩（式26-27, Delta法同权）
%  ========================================================================
fprintf('\n[节点功率] 计算节点注入有功功率统计矩...\n');

% MC基准: 直接计算 P_i = Re(V_i * conj(sum_j Y_ij * V_j))
V_mc_mat = V_mc_all.';                          % [nc x N] complex
I_mc_mat = Y3ph_ordered * V_mc_mat;              % [nc x N] current injection
Node_Power_P_direct = real(V_mc_mat .* conj(I_mc_mat));  % [nc x N]
mean_NodeP_mc = mean(Node_Power_P_direct, 2);
var_NodeP_mc  = var(Node_Power_P_direct, 0, 2);

% 二次型MC验证: tr(M_i * C) 应与直接法均值一致, C = (1/N)*V_real*V_real'
V_mc_T = v_vec';  % [2*nc x N]
C_data = (1/N) * (V_mc_T * V_mc_T');            % [2nc x 2nc] 二阶矩矩阵
NodeP_quad_check = zeros(nc, 1);
for i = 1:nc
    NodeP_quad_check(i) = sum(Node_P_matrix{i}(:) .* C_data(:));
end
thr_quad = 0.01 * max(abs(mean_NodeP_mc));
quad_vs_direct = max(abs(NodeP_quad_check - mean_NodeP_mc) ./ ...
    max(abs(mean_NodeP_mc), thr_quad)) * 100;
fprintf('[节点功率] 二次型-直接法均值一致性: %.4f%%\n', quad_vs_direct);

% 解析计算（式26-27）
tic_nodeP = tic;
mean_NodeP_com = zeros(nc, gmm_num);
var_NodeP_com  = zeros(nc, gmm_num);
for i = 1:nc
    for k = 1:gmm_num
        Sigma_k = gmm_fit.Sigma(:,:,k);
        mu_k = gmm_fit.mu(k,:)';
        mean_NodeP_com(i,k) = trace(Node_P_matrix{i} * Sigma_k) ...
            + mu_k' * Node_P_matrix{i} * mu_k;
        var_NodeP_com(i,k) = 2*trace((Node_P_matrix{i}*Sigma_k)^2) ...
            + 4*mu_k'*Node_P_matrix{i}*Sigma_k*Node_P_matrix{i}*mu_k;
    end
end
mean_NodeP_this = mean_NodeP_com * com_pi';
var_NodeP_this  = var_NodeP_com * com_pi' + ...
    sum(com_pi .* (mean_NodeP_com - mean_NodeP_this).^2, 2);
t_nodeP_stat = toc(tic_nodeP);

% 误差（分母保护: 与电压同口径，用 max 的 1% 兜底）
sig_thr = 0.01 * max(abs(mean_NodeP_mc));
sig_idx = abs(mean_NodeP_mc) > sig_thr;
thr_NP_m = sig_thr;
thr_NP_v = 0.01 * max(abs(var_NodeP_mc));
error_NodeP_mean = abs(mean_NodeP_mc - mean_NodeP_this) ./ max(abs(mean_NodeP_mc), thr_NP_m) * 100;
error_NodeP_var  = abs(var_NodeP_mc - var_NodeP_this)   ./ max(abs(var_NodeP_mc),  thr_NP_v) * 100;
error_NodeP_mean(~isfinite(error_NodeP_mean)) = 0;
error_NodeP_var(~isfinite(error_NodeP_var))   = 0;

fprintf('[节点功率] 显著注入节点: %d/%d\n', sum(sig_idx), nc);
fprintf('[节点功率] P均值最大误差(显著节点): %.4f%%, P方差最大误差(显著节点): %.4f%%\n', ...
    max(error_NodeP_mean(sig_idx)), max(error_NodeP_var(sig_idx)));
fprintf('[节点功率] P均值最大误差(全网): %.4f%%, P方差最大误差(全网): %.4f%%\n', ...
    max(error_NodeP_mean), max(error_NodeP_var));
fprintf('[节点功率] 统计矩计算耗时: %.4f s\n', t_nodeP_stat);

%% ========================================================================
%  段落10: PEM (Hong 2m+1) 对比方法
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  段落10: PEM (Hong 2m+1) — IEEE 123\n');
fprintf('%s\n', repmat('=',1,70));

m_rv = 24;
pemCacheFile = fullfile(cacheDir, sprintf('PEM_3ph_%s_T%d.mat', caseName, N));

if isfile(pemCacheFile)
    fprintf('[PEM] 发现缓存: %s\n', pemCacheFile);
    pemData = load(pemCacheFile);
    V_pem_all    = pemData.V_pem_all;
    w_pem        = pemData.w_pem;
    w0_pem       = pemData.w0_pem;
    V_pem_center = pemData.V_pem_center;
    t_pem_dss    = pemData.t_pem_dss;
    mu_X         = pemData.mu_X;
    sigma_X      = pemData.sigma_X;
    fprintf('[PEM] 缓存加载完成\n');
else
    fprintf('[PEM] 无缓存，执行 2m+1=%d 次OpenDSS潮流...\n', 2*m_rv+1);

    X_data = zeros(N, m_rv);
    for c = 1:m_rv
        X_data(:, c) = load_data{:, c};
    end
    mu_X    = mean(X_data);
    sigma_X = std(X_data);
    skew_X  = skewness(X_data);

    % Hong 2m+1 concentration points & weights
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

    DSSObj = actxserver('OpenDSSEngine.DSS');
    assert(DSSObj.Start(0), '[OpenDSS] PEM启动失败');

    tic;
    V_pem_center = tools.eval_opendss_single(DSSObj, dssFile, mu_X(:), ...
        load_sequence(:), loadKW_base, loadKvar_base, nLoad);

    V_pem_all = zeros(length(V_pem_center), 2*m_rv);
    w_pem = zeros(2*m_rv, 1);
    h_pem = waitbar(0, 'PEM: 0/48', 'Name', 'IEEE123 PEM 2m+1');
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

    DSSObj.delete; clear DSSObj;

    save(pemCacheFile, 'V_pem_all', 'w_pem', 'w0_pem', 'V_pem_center', ...
        't_pem_dss', 'mu_X', 'sigma_X', '-v7.3');
    fprintf('[PEM] 缓存已保存: %s\n', pemCacheFile);
end

% PEM voltage magnitude mean/variance
Vmag_pem_all = abs(V_pem_all);          % [nc_pem x 2m]
Vmag_pem_center = abs(V_pem_center);    % [nc_pem x 1]
nc_pem = length(V_pem_center);

mean_Vmag_pem = w0_pem * Vmag_pem_center + Vmag_pem_all * w_pem;
var_Vmag_pem  = w0_pem * Vmag_pem_center.^2 + Vmag_pem_all.^2 * w_pem ...
              - mean_Vmag_pem.^2;

% PEM errors vs MCS (reuse mean_Vmag_mc, var_Vmag_mc from 段落8)
error_Vmag_mean_pem = abs(mean_Vmag_pem(1:nc) - mean_Vmag_mc) ...
    ./ abs(mean_Vmag_mc) * 100;
error_Vmag_var_pem  = abs(var_Vmag_pem(1:nc) - var_Vmag_mc) ...
    ./ abs(var_Vmag_mc) * 100;
error_Vmag_mean_pem(~isfinite(error_Vmag_mean_pem)) = 0;
error_Vmag_var_pem(~isfinite(error_Vmag_var_pem))   = 0;

fprintf('[PEM] Vmag均值最大误差: %.4f%%, 方差最大误差: %.4f%%\n', ...
    max(error_Vmag_mean_pem), max(error_Vmag_var_pem));
fprintf('[PEM] 耗时: %.2f s\n', t_pem_dss);

%% ========================================================================
%  目标产物: Fig.7, Fig.7b, Fig.8, Tab.3
%  ========================================================================
fprintf('\n[产物] 生成IEEE 123目标图表...\n');

% === Fig.7: 全网电压幅值误差散点图 ===
fprintf('\n[图7] 全网电压幅值误差散点图...\n');
figure('Units','pixels','Position',[100,100,1000,500]);
hold on;
scatter(1:nc, error_Vmag_mean, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(1:nc, error_Vmag_var,  30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
hold off;
set(gca, 'fontsize', fontSZ, 'fontname', fontName_en, 'box', 'off');
xlabel('节点编号', 'FontName', fontName_cn, 'FontSize', fontSZ);
ylabel('误差 (%)', 'FontName', fontName_cn, 'FontSize', fontSZ);
legend({'均值误差','方差误差'}, 'FontName', fontName_cn, 'FontSize', fontSZ, 'Location', 'best');
title(sprintf('IEEE 123 全网电压幅值误差 (K=%d)', gmm_num), ...
    'FontName', fontName_cn, 'FontSize', fontSZ_title);
saveas(gcf, fullfile(resultDir, 'Fig07_IEEE123_VmagError_Scatter.png'));
fprintf('[图7] 已保存\n');

% === Fig.7b: 全网有功功率误差散点图 ===
fprintf('\n[图7b] 全网有功功率误差散点图...\n');
figure('Units','pixels','Position',[150,150,1000,500]);
hold on;
scatter(1:nc, error_NodeP_mean, 30, 'b', 'filled', 'MarkerFaceAlpha', 0.6);
scatter(1:nc, error_NodeP_var,  30, 'r', 'filled', 'MarkerFaceAlpha', 0.6);
hold off;
set(gca, 'fontsize', fontSZ, 'fontname', fontName_en, 'box', 'off');
xlabel('节点编号', 'FontName', fontName_cn, 'FontSize', fontSZ);
ylabel('误差 (%)', 'FontName', fontName_cn, 'FontSize', fontSZ);
legend({'均值误差','方差误差'}, 'FontName', fontName_cn, 'FontSize', fontSZ, 'Location', 'best');
title(sprintf('IEEE 123 全网有功功率误差 (K=%d)', gmm_num), ...
    'FontName', fontName_cn, 'FontSize', fontSZ_title);
saveas(gcf, fullfile(resultDir, 'Fig07b_IEEE123_PowerError_Scatter.png'));
fprintf('[图7b] 已保存\n');

% === Fig.8: K值敏感性双Y轴图 ===
fprintf('\n[图8] K值敏感性分析...\n');
K_values = [5, 10, 15, 20, 30, 40, 50, 60, 75, 100];
max_var_err = zeros(length(K_values), 1);
comp_time_K = zeros(length(K_values), 1);

n_gncs_samples_K = 1e4;

for ki = 1:length(K_values)
    K_cur = K_values(ki);
    fprintf('  K=%d (%d/%d)...', K_cur, ki, length(K_values));

    tic_k = tic;

    try
        gmm_k = fitgmdist(v_vec, K_cur, 'RegularizationValue', 1e-6, ...
            'Options', statset('MaxIter', 300, 'TolFun', 1e-6), ...
            'CovarianceType', 'full', 'SharedCovariance', false);
    catch ME
        fprintf(' GMM拟合失败: %s\n', ME.message);
        max_var_err(ki) = NaN;
        comp_time_K(ki) = NaN;
        continue;
    end

    pi_k = gmm_k.ComponentProportion;
    B_k = cell(K_cur, 1);
    for kk = 1:K_cur
        Sk = gmm_k.Sigma(:,:,kk);
        [Vk_eig, Dk_eig] = eig(Sk);
        dk = diag(Dk_eig); dk(dk<0) = 0;
        B_k{kk} = Vk_eig * diag(sqrt(dk));
    end

    mean_Vmag_k = zeros(nc, 1);
    var_Vmag_k  = zeros(nc, 1);
    for i = 1:nc
        smp = [];
        for kk = 1:K_cur
            BMB = B_k{kk}' * M_Vmag2{i} * B_k{kk};
            BMB_sym = 0.5*(BMB + BMB');
            [Qv, Dv] = eig(BMB_sym);
            lv = diag(Dv);
            mu_kk = gmm_k.mu(kk,:)';
            av = 2*B_k{kk}'*M_Vmag2{i}*mu_kk;
            bv = Qv'*av;
            alpha_v = mu_kk'*M_Vmag2{i}*mu_kk;
            tol = 1e-6*max(abs(lv));
            inz = abs(lv) > tol;
            cv = zeros(size(lv)); cv(inz) = bv(inz)./(2*lv(inz));
            dv = cv.^2;
            lv_nz = lv(inz); dv_nz = dv(inz);
            c_nz = alpha_v - sum(lv_nz.*dv_nz);

            nk = round(n_gncs_samples_K * pi_k(kk));
            if nk == 0, continue; end
            Z = randn(length(lv_nz), nk);
            sq = c_nz + sum(lv_nz.*(Z + sqrt(dv_nz)).^2, 1);
            smp = [smp, sqrt(max(sq, 0))];
        end
        if ~isempty(smp)
            mean_Vmag_k(i) = mean(smp);
            var_Vmag_k(i)  = var(smp);
        end
    end

    comp_time_K(ki) = toc(tic_k);
    err_var_k = abs((var_Vmag_mc - var_Vmag_k) ./ max(abs(var_Vmag_mc), thr_vmag_var)) * 100;
    max_var_err(ki) = max(err_var_k);

    fprintf(' maxVarErr=%.2f%%, time=%.1fs\n', max_var_err(ki), comp_time_K(ki));
end

valid = ~isnan(max_var_err);
figure('Units','pixels','Position',[100,100,900,500]);
yyaxis left;
bar(find(valid), max_var_err(valid), 0.5, 'FaceColor', [0 0.45 0.74], 'FaceAlpha', 0.7);
ylabel('最大方差误差 (%)', 'FontName', fontName_cn, 'FontSize', fontSZ);
set(gca, 'XTick', 1:length(K_values), 'XTickLabel', arrayfun(@num2str, K_values, 'UniformOutput', false));

yyaxis right;
plot(find(valid), comp_time_K(valid), '-o', 'Color', [0.85 0.33 0.1], 'LineWidth', 2, 'MarkerSize', 8);
ylabel('计算时间 (s)', 'FontName', fontName_cn, 'FontSize', fontSZ);

xlabel('GMM分量数 K', 'FontName', fontName_cn, 'FontSize', fontSZ);
set(gca, 'fontsize', fontSZ, 'fontname', fontName_en, 'box', 'off');
title('IEEE 123 K值敏感性分析', 'FontName', fontName_cn, 'FontSize', fontSZ_title);
legend({'方差误差','计算时间'}, 'FontName', fontName_cn, 'FontSize', fontSZ, 'Location', 'best');
saveas(gcf, fullfile(resultDir, 'Fig08_IEEE123_KSensitivity.png'));
fprintf('[图8] 已保存\n');

% === 表3: 计算时间拆解与对比 ===
t_common = t_gmm + t_gncs + t_vmag_gncs;
t_total_exact = t_common + t_vmag_exact + t_stat + t_nodeP_stat;
t_total_delta = t_common + t_vmag_delta + t_stat + t_nodeP_stat;

fprintf('\n');
fprintf('================================================================================\n');
fprintf('  表3: IEEE 123 概率潮流计算时间拆解与对比 (K=%d, N=%d)\n', gmm_num, N);
fprintf('================================================================================\n');
fprintf('  %-32s | %10s | %s\n', '方法 / 阶段', '耗时 (s)', '备注');
fprintf('  ---------------------------------|------------|---------------------------\n');
fprintf('  %-32s | %10.2f | %s\n', '本文-GMM拟合',              t_gmm,       '一次性');
fprintf('  %-32s | %10.2f | %s\n', '本文-二次型参数构建(功率)',   t_gncs,       '一次性');
fprintf('  %-32s | %10.2f | %s\n', '本文-电压幅值GNCS参数',      t_vmag_gncs, '一次性');
fprintf('  %-32s | %10.2f | %s\n', '本文-功率统计矩(解析)',       t_stat + t_nodeP_stat, '式26-27');
fprintf('  %-32s | %10.2f | %s\n', '本文-电压统计矩(精确Laplace)', t_vmag_exact, 'GNCS采样');
fprintf('  %-32s | %10.4f | %s\n', '本文-电压统计矩(Delta近似)',   t_vmag_delta, '式42-46b');
fprintf('  ---------------------------------|------------|---------------------------\n');
fprintf('  %-32s | %10.2f | %s\n', '本文总耗时(精确Laplace)',    t_total_exact, '');
fprintf('  %-32s | %10.2f | %s\n', '本文总耗时(Delta近似)',      t_total_delta, '');
fprintf('  ---------------------------------|------------|---------------------------\n');
if t_mc > 0
    fprintf('  %-32s | %10.2f | %d次潮流\n', 'MCS总耗时', t_mc, N);
    fprintf('  %-32s | %10.1f | \n', '加速比(MCS/Delta)', t_mc / t_total_delta);
end
fprintf('  %-32s | %10.2f | %d次潮流\n', 'PEM总耗时', t_pem_dss, 2*m_rv+1);
fprintf('================================================================================\n');

%% ========================================================================
%  汇总
%  ========================================================================
fprintf('\n%s\n', repmat('=',1,70));
fprintf('  AQPPF_3ph_123bus 计算完成\n');
fprintf('%s\n', repmat('=',1,70));
fprintf('  算例: %s, 节点数: %d, K=%d\n', caseName, nc, gmm_num);
fprintf('  电压幅值均值最大误差: 精确=%.4f%%, Delta=%.4f%%, PEM=%.4f%%\n', ...
    max(error_Vmag_mean), max(error_Vmag_mean_delta), max(error_Vmag_mean_pem));
fprintf('  电压幅值方差最大误差: 精确=%.4f%%, Delta=%.4f%%, PEM=%.4f%%\n', ...
    max(error_Vmag_var), max(error_Vmag_var_delta), max(error_Vmag_var_pem));
fprintf('  计算时间: 精确=%.1fs, Delta=%.1fs, PEM=%.1fs\n', ...
    t_total_exact, t_total_delta, t_pem_dss);
if t_mc > 0
    fprintf('  计算时间: MCS=%.1fs\n', t_mc);
end
fprintf('  结果目录: %s\n', resultDir);
fprintf('%s\n', repmat('=',1,70));

wsFile = fullfile(resultDir, 'workspace_123bus.mat');
if isfile(wsFile), delete(wsFile); end
save(wsFile, '-v7.3');
fprintf('[保存] 工作空间已保存\n');
