function [Mbus_Re,Mbus_Im] = get_node_power_matrix(mpc)
% get_node_power_matrix 计算潮流注入功率二次型参数
% [Mbus_Re,Mbus_Im] = get_node_power_matrix(mpc)

% 转换为内部编号，防止case18这样的奇葩case
mpc = ext2int(mpc);
% 计算导纳矩阵
Ybus = makeYbus(mpc);
Mbus_Re = {}; % 用于计算节点有功的系数
Mbus_Im = {}; % 用于计算节点无功的系数
for i = 1:length(Ybus)
    ReA = zeros(length(Ybus)*2); % 对应有功功率
    ImA = zeros(length(Ybus)*2); % 对应无功功率

    Y_e = real(Ybus(:,i));
    Y_f = imag(Ybus(:,i));

    % P_i = Σ_j [G_ij*(e_i*e_j+f_i*f_j) + B_ij*(f_i*e_j-e_i*f_j)]
    % M_P = [G_i, -B_i; B_i, G_i]
    ReA(i,:) = [Y_e;-Y_f]; 
    ReA(i+length(Ybus),:) = [Y_f;Y_e];
    Mbus_Re{i} = ReA;  % 用元胞数组来储每个节点对应的系数矩阵

    % Q_i = Σ_j [G_ij*(f_i*e_j-e_i*f_j) - B_ij*(e_i*e_j+f_i*f_j)]
    % M_Q = [-B_i, -G_i; G_i, -B_i]
    ImA(i,:) = [-Y_f;-Y_e];
    ImA(i+length(Ybus),:) = [Y_e;-Y_f];
    Mbus_Im{i} = ImA;
end
end

