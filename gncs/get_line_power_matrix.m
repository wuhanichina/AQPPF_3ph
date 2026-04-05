function [Aline_Re, Aline_Im] = get_line_power_matrix(mpc)

Aline_Re = {};% 用于计算线路有功的系数
Aline_Im = {};% 用于计算线路无功的系数
Aloss_Re = {};% 用于计算线路网损的系数

% 转换为内部编号，防止case18这样的奇葩case
mpc = ext2int(mpc);
% 计算导纳矩阵
Ybus = makeYbus(mpc);
nl = length(mpc.branch(mpc.branch(:, 11)==1, 1)); % 线路数量
nb = size(Ybus, 1); % 节点数量
    
for l = 1:nl % 对每条线路执行操作
    % 获取线路两端的节点编号
    i = mpc.branch(l,1);
    j = mpc.branch(l,2);

    % Ysc = 1 ./ (mpc.branch(1:32, 3) + 1j * mpc.branch(1:32, 4));
    % 获取线路导纳
    G = real(Ybus(i, j));
    B = imag(Ybus(i, j));

    % 初始化系数矩阵
    M_P = zeros(2*nb, 2*nb); % 有功功率
    M_Q = zeros(2*nb, 2*nb); % 无功功率

    % 线路原始导纳 y_ij = -Y_bus(i,j)
    % P_ij = -G*(|Vi|^2 - Re(Vi*Vj')) + B*Im(Vi*Vj')
    % Q_ij = B*(|Vi|^2 - Re(Vi*Vj')) + G*Im(Vi*Vj')

    % 填充有功功率的二次型系数矩阵
    M_P(i, i) = -G;
    M_P(i, j) = G;
    M_P(nb+i, nb+i) = -G;
    M_P(nb+i, nb+j) = G;
    M_P(i, nb+j) = -B;
    M_P(nb+i, j) = B;

    % 填充无功功率的二次型系数矩阵
    M_Q(i, i) = B;
    M_Q(i, j) = -B;
    M_Q(nb+i, nb+i) = B;
    M_Q(nb+i, nb+j) = -B;
    M_Q(i, nb+j) = -G;
    M_Q(nb+i, j) = G;

    Aline_Re{l} = M_P; % 用元胞数组来存储每个节点对应的系数矩阵
    Aline_Im{l} = M_Q; % 用元胞数组来存储每个节点对应的系数矩阵
end
end