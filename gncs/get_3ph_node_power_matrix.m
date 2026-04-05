function [M_P, M_Q] = get_3ph_node_power_matrix(Y3ph, nc)
% GET_3PH_NODE_POWER_MATRIX 构建三相不平衡配电网节点注入功率的二次型参数矩阵
%
% 语法:
%   [M_P, M_Q] = get_3ph_node_power_matrix(Y3ph, nc)
%
% 输入:
%   Y3ph - 三相节点导纳矩阵 [nc x nc] (复数)
%          nc 为计算节点总数（各相分别计数）
%   nc   - 计算节点总数
%
% 输出:
%   M_P  - cell{nc,1}，M_P{i} 为节点i注入有功功率的二次型参数矩阵 [2*nc x 2*nc]
%   M_Q  - cell{nc,1}，M_Q{i} 为节点i注入无功功率的二次型参数矩阵 [2*nc x 2*nc]
%
% 说明:
%   对于实同构后的电压向量 V = [V_e; V_f] (2*nc x 1),
%   节点i的注入功率表示为:
%     P_i = V' * M_P{i} * V
%     Q_i = V' * M_Q{i} * V
%
%   其中 M_P{i} 和 M_Q{i} 的结构为（主要内容.md 式5）:
%     M_P{i} = [G_i, -B_i; B_i, G_i]
%     M_Q{i} = [-B_i, -G_i; G_i, -B_i]
%
%   子矩阵 G_i 和 B_i 从三相导纳矩阵的第i行提取:
%     G_i(i,j) = Re(Y3ph(i,j)),  B_i(i,j) = Im(Y3ph(i,j))
%     （仅第i行非零，其余行全零）
%
%   最终对称化处理（主要内容.md 式6）:
%     M_P{i} = 0.5 * (M_P{i} + M_P{i}')
%     M_Q{i} = 0.5 * (M_Q{i} + M_Q{i}')
%
% 版本: v1.1 (修正: 提取第i行而非第i列，兼容非对称Y)
% 日期: 2026-02-26

    arguments
        Y3ph (:,:) {mustBeNonempty}
        nc   (1,1) double {mustBePositive, mustBeInteger}
    end

    assert(size(Y3ph,1) == nc && size(Y3ph,2) == nc, ...
        '导纳矩阵维度 [%d x %d] 与计算节点数 nc=%d 不匹配', ...
        size(Y3ph,1), size(Y3ph,2), nc);

    M_P = cell(nc, 1);
    M_Q = cell(nc, 1);

    for i = 1:nc
        % 提取导纳矩阵第i行的实部和虚部
        Y_row = Y3ph(i, :);              % 第i行 [1 x nc]
        Y_e = real(Y_row);               % G_{i,:} [1 x nc]
        Y_f = imag(Y_row);               % B_{i,:} [1 x nc]

        % 构建有功功率二次型矩阵（主要内容.md 式5）
        %   P_i = Σ_j [G_ij*(e_i*e_j + f_i*f_j) + B_ij*(f_i*e_j - e_i*f_j)]
        %   M_P{i} = [G_i, -B_i; B_i, G_i]
        % G_i 和 B_i 均为 nc x nc 稀疏矩阵，仅第 i 行非零
        ReA = zeros(2*nc, 2*nc);
        ReA(i, 1:nc)         = Y_e;      % G_row_i
        ReA(i, nc+1:2*nc)    = -Y_f;     % -B_row_i
        ReA(nc+i, 1:nc)      = Y_f;      % B_row_i
        ReA(nc+i, nc+1:2*nc) = Y_e;      % G_row_i

        % 构建无功功率二次型矩阵（主要内容.md 式5）
        %   Q_i = Σ_j [G_ij*(f_i*e_j - e_i*f_j) - B_ij*(e_i*e_j + f_i*f_j)]
        %   M_Q{i} = [-B_i, -G_i; G_i, -B_i]
        ImA = zeros(2*nc, 2*nc);
        ImA(i, 1:nc)         = -Y_f;     % -B_row_i
        ImA(i, nc+1:2*nc)    = -Y_e;     % -G_row_i
        ImA(nc+i, 1:nc)      = Y_e;      % G_row_i
        ImA(nc+i, nc+1:2*nc) = -Y_f;     % -B_row_i

        % 对称化处理（主要内容.md 式6）
        M_P{i} = 0.5 * (ReA + ReA');
        M_Q{i} = 0.5 * (ImA + ImA');
    end

end
