function M_I = get_line_current(mpc)

% 转换为内部编号，防止case18这样的奇葩case
mpc = ext2int(mpc);
% 计算导纳矩阵
Ybus = makeYbus(mpc);
nb = size(Ybus, 1); % 节点数量

M_I = zeros(2*nb, 2*nb);
% Fill the matrix for current injections
for i = 1:nb
    for j = 1:nb
        M_I(i, j) = real(Ybus(i, j));                % Real part affects real current
        M_I(i, nb + j) = -imag(Ybus(i, j));          % Imaginary part affects real current
        M_I(nb + i, j) = imag(Ybus(i, j));           % Real part affects imaginary current
        M_I(nb + i, nb + j) = real(Ybus(i, j));      % Imaginary part affects imaginary current
    end
end

end