function [M_line_P, M_line_Q] = get_3ph_line_power_matrix(Y3ph, branchData, nc)
% GET_3PH_LINE_POWER_MATRIX 构建三相配电网线路首端功率的二次型参数矩阵
%
% 语法:
%   [M_line_P, M_line_Q] = get_3ph_line_power_matrix(Y3ph, branchData, nc)
%
% 输入:
%   Y3ph       - 三相节点导纳矩阵 [nc x nc] (复数)
%   branchData - 线路连接信息 [nLine x 2], 每行为 [from_node, to_node]
%                使用计算节点编号 (1-based)
%   nc         - 计算节点总数
%
% 输出:
%   M_line_P - cell{nLine,1}, 各线路首端有功功率二次型矩阵 [2*nc x 2*nc]
%   M_line_Q - cell{nLine,1}, 各线路首端无功功率二次型矩阵 [2*nc x 2*nc]
%
% 说明:
%   线路l连接节点i到j，其首端功率利用线路原始导纳 y_ij = -Y_bus(i,j):
%     P_ij = Re(V_i * conj(y_ij*(V_i - V_j)))
%          = -G*(|Vi|^2 - Re(Vi*Vj')) + B*Im(Vi*Vj')
%     Q_ij = Im(V_i * conj(y_ij*(V_i - V_j)))
%          = B*(|Vi|^2 - Re(Vi*Vj')) + G*Im(Vi*Vj')
%   其中 G = Re(Y_bus(i,j)), B = Im(Y_bus(i,j))
%
%   对于三相线路，各相功率还包含互耦项的贡献，这些通过
%   get_3ph_line_power_matrix_coupled 函数处理。
%
% 版本: v2.0 (修正符号错误)
% 日期: 2026-02-08

    arguments
        Y3ph       (:,:) {mustBeNonempty}
        branchData (:,2) double {mustBePositive, mustBeInteger}
        nc         (1,1) double {mustBePositive, mustBeInteger}
    end

    nLine = size(branchData, 1);
    M_line_P = cell(nLine, 1);
    M_line_Q = cell(nLine, 1);

    for l = 1:nLine
        i = branchData(l, 1);  % 首端节点
        j = branchData(l, 2);  % 末端节点

        % Y_bus 导纳矩阵元素
        G = real(Y3ph(i, j));
        B = imag(Y3ph(i, j));

        % 线路原始导纳 y_ij = -Y_bus(i,j) = (-G) + j(-B)
        % P_ij = -G*(|Vi|^2 - Re(Vi*Vj')) + B*Im(Vi*Vj')  (自导纳项)

        % 有功功率二次型矩阵
        M_P = zeros(2*nc, 2*nc);
        M_P(i, i)       = -G;
        M_P(i, j)       = G;
        M_P(nc+i, nc+i) = -G;
        M_P(nc+i, nc+j) = G;
        M_P(i, nc+j)    = -B;
        M_P(nc+i, j)    = B;

        % 无功功率二次型矩阵
        % Q_ij = B*(|Vi|^2 - Re(Vi*Vj')) + G*Im(Vi*Vj')
        M_Q = zeros(2*nc, 2*nc);
        M_Q(i, i)       = B;
        M_Q(i, j)       = -B;
        M_Q(nc+i, nc+i) = B;
        M_Q(nc+i, nc+j) = -B;
        M_Q(i, nc+j)    = -G;
        M_Q(nc+i, j)    = G;

        % 对称化
        M_line_P{l} = 0.5 * (M_P + M_P');
        M_line_Q{l} = 0.5 * (M_Q + M_Q');
    end

end
