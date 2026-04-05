function [M_loss_P, M_loss_Q] = get_system_loss_matrices(mpc)
% 计算全网有功和无功网损的二次型系数矩阵
% Ybus - 系统导纳矩阵
% n - 节点总数
% 转换为内部编号，防止case18这样的奇葩case
mpc = ext2int(mpc);
% 计算导纳矩阵
Ybus = makeYbus(mpc);
n = size(Ybus, 1);
% 初始化全网网损系数矩阵
M_loss_P = zeros(2*n, 2*n);
M_loss_Q = zeros(2*n, 2*n);

for i = 1:n
    for j = i+1:n
        if i ~= j && Ybus(i, j) ~= 0
            Y = Ybus(i, j);
            G = real(Y);
            B = imag(Y);

            % 有功网损系数
            M_loss_P(i, i) = M_loss_P(i, i) + G;
            M_loss_P(i, j) = M_loss_P(i, j) - G;
            M_loss_P(j, i) = M_loss_P(j, i) - G;
            M_loss_P(j, j) = M_loss_P(j, j) + G;

            M_loss_P(n+i, n+i) = M_loss_P(n+i, n+i) + G;
            M_loss_P(n+i, n+j) = M_loss_P(n+i, n+j) - G;
            M_loss_P(n+j, n+i) = M_loss_P(n+j, n+i) - G;
            M_loss_P(n+j, n+j) = M_loss_P(n+j, n+j) + G;

            % 无功网损系数
            M_loss_Q(i, i) = M_loss_Q(i, i) - B;
            M_loss_Q(i, j) = M_loss_Q(i, j) + B;
            M_loss_Q(j, i) = M_loss_Q(j, i) + B;
            M_loss_Q(j, j) = M_loss_Q(j, j) - B;

            M_loss_Q(n+i, n+i) = M_loss_Q(n+i, n+i) - B;
            M_loss_Q(n+i, n+j) = M_loss_Q(n+i, n+j) + B;
            M_loss_Q(n+j, n+i) = M_loss_Q(n+j, n+i) + B;
            M_loss_Q(n+j, n+j) = M_loss_Q(n+j, n+j) - B;
        end
    end
end
M_loss_P = -M_loss_P;
M_loss_Q = -M_loss_Q;
end