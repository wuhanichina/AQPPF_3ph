function Aloss_Re = get_line_loss_matrix(mpc)

% 用于计算线路网损的系数
Aloss_Re = {};
% 转换为内部编号，防止case18这样的奇葩case
mpc = ext2int(mpc);
% 计算导纳矩阵
Ybus = makeYbus(mpc);

% Get the number of buses
nb = size(Ybus, 1);
% Get the number of lines
nl = length(mpc.branch(mpc.branch(:, 11)==1, 1)); 

for l=1:nl
% Get nodes connected by this line
from = mpc.branch(l, 1);
to = mpc.branch(l, 2);
G = real(Ybus(from, to));  % Get the conductance between these nodes

% Initialize the loss matrix for this specific line
M_loss = zeros(2*nb, 2*nb);

% Update the loss matrix for the real parts
M_loss(from, from) = G;
M_loss(from, to) = -G;
M_loss(to, from) = -G;
M_loss(to, to) = G;

% Update the loss matrix for the imaginary parts
M_loss(nb + from, nb + from) = G;
M_loss(nb + from, nb + to) = -G;
M_loss(nb + to, nb + from) = -G;
M_loss(nb + to, nb + to) = G;

Aloss_Re{l} = -M_loss;
end

end