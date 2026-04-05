function M_loss = get_3ph_line_loss_matrix(Y3ph, branchData, nc)
% GET_3PH_LINE_LOSS_MATRIX 构建三相配电网线路有功损耗的二次型参数矩阵
%
% 语法:
%   M_loss = get_3ph_line_loss_matrix(Y3ph, branchData, nc)
%
% 输入:
%   Y3ph       - 三相节点导纳矩阵 [nc x nc] (复数)
%   branchData - 线路连接信息 [nLine x 2], 每行为 [from_node, to_node]
%   nc         - 计算节点总数
%
% 输出:
%   M_loss - cell{nLine,1}, 各线路有功损耗二次型矩阵 [2*nc x 2*nc]
%
% 说明:
%   线路l的有功损耗 = P_loss = V' * M_loss{l} * V
%   利用线路两端的电导 G 构建:
%     M_loss 在实部和虚部块中均为 G*([i,i]-[i,j]-[j,i]+[j,j]) 形式
%
% 版本: v1.0
% 日期: 2026-02-07

    arguments
        Y3ph       (:,:) {mustBeNonempty}
        branchData (:,2) double {mustBePositive, mustBeInteger}
        nc         (1,1) double {mustBePositive, mustBeInteger}
    end

    nLine = size(branchData, 1);
    M_loss = cell(nLine, 1);

    for l = 1:nLine
        from = branchData(l, 1);
        to   = branchData(l, 2);

        G = real(Y3ph(from, to));

        M_l = zeros(2*nc, 2*nc);

        % 实部块 (上左 nc x nc)
        M_l(from, from)  = G;
        M_l(from, to)    = -G;
        M_l(to, from)    = -G;
        M_l(to, to)      = G;

        % 虚部块 (下右 nc x nc)
        M_l(nc+from, nc+from) = G;
        M_l(nc+from, nc+to)   = -G;
        M_l(nc+to, nc+from)   = -G;
        M_l(nc+to, nc+to)     = G;

        % 取负号 (与AQPPF中get_line_loss_matrix一致)
        M_loss{l} = -M_l;
    end

end
