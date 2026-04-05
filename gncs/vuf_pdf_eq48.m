function f = vuf_pdf_eq48(y, lam_pos, lam_neg, rho, sig_ratio)
%VUF_PDF_EQ48 VUF概率密度函数（主要内容.md 式47-51，带大λ+近似）
%
% 根据主要内容.md 附录A.2推导，VUF的PDF基于两个相关非中心χ²(df=2)变量之比的分布。
%   R = |V-|²/|V+|² = (σ²-/σ²+) · (W-/W+)
%   其中 W± ~ χ²(2, λ±) 为标准化的非中心卡方变量
%   VUF = 100√R (%)
%
% 当 λ+ 很大时（实际电力系统中 |V+| >> σ+），直接计算式(48)会因指数下溢
% 而得到零。此时使用Rice近似：将 |V+| 视为近似常量，
% VUF ≈ 100·sig_ratio/√λ+ · |V-|/σ-，其中 |V-|/σ- 服从Rice分布。
%
% 输入:
%   y         - VUF值向量 (%), 任意形状
%   lam_pos   - 标准化正序非中心参数 λ+ = |E[V+]|²/σ²+ (scalar)
%   lam_neg   - 标准化负序非中心参数 λ- = |E[V-]|²/σ²- (scalar)
%   rho       - |V+|² 与 |V-|² 的相关系数 (scalar, |ρ| < 1)
%   sig_ratio - σ-/σ+ 尺度比 (scalar, default=1)
%
% 输出:
%   f - VUF的概率密度值（与y同维）

    arguments
        y         double
        lam_pos   (1,1) double {mustBeNonnegative}
        lam_neg   (1,1) double {mustBeNonnegative}
        rho       (1,1) double
        sig_ratio (1,1) double {mustBePositive} = 1.0
    end

    sz = size(y);
    y = y(:);  % 统一为列向量

    f = zeros(size(y));
    valid = y > 0;
    if ~any(valid)
        f = reshape(f, sz);
        return;
    end
    yv = y(valid);

    % 判断是否需要使用大 λ+ 近似
    % 当 λ+/(2(1+|ρ|)) > 500 时，exp(-)项必定下溢为零
    use_rice = (lam_pos / (2*(1+abs(rho))) > 500);

    if use_rice
        % ===== Rice近似（大λ+渐近）=====
        % 当 λ+ >> 1 时, W+ ≈ λ+ + 2 (近似常量)
        % VUF = 100·sig_ratio·√(W-/W+) ≈ 100·sig_ratio/√(λ++2)·√W-
        % 令 c = 100·sig_ratio/√(λ++2), t = y/c = √W-
        % √W- 的PDF为 Rice: f_t(t) = t·exp(-(t²+λ-)/2)·I₀(t√λ-)
        % f_VUF(y) = (1/c)·f_t(y/c)
        c = 100 * sig_ratio / sqrt(lam_pos + 2);
        t = yv / c;
        sqrt_lam_neg = sqrt(lam_neg);
        bess_arg = t * sqrt_lam_neg;

        % 数值稳定: exp(-(t²+λ-)/2)·I₀(bess_arg) =
        %   exp(-(t²+λ-)/2 + bess_arg)·besseli(0, bess_arg, 1)
        log_arg = -(t.^2 + lam_neg)/2 + bess_arg;
        bess_scaled = besseli(0, bess_arg, 1);

        f_rice = (t / c) .* exp(log_arg) .* bess_scaled;
        f(valid) = max(f_rice, 0);
    else
        % ===== 原始式(48) =====
        one_minus_rho2 = 1 - rho^2;
        assert(one_minus_rho2 > 0, 'VUF PDF: |ρ| 必须 < 1, 当前 ρ = %.6f', rho);

        y_std = yv / sig_ratio;
        r = (y_std / 100).^2;
        sqrt_r = y_std / 100;
        sqrt_lam_prod = sqrt(lam_pos * lam_neg);

        denom = 2 * one_minus_rho2 * (1 + r);

        exp_num = -(lam_pos + r * lam_neg - 2 * rho * sqrt_r * sqrt_lam_prod);
        exp_arg = exp_num ./ denom;

        bess_num = rho * (lam_pos + r * lam_neg) - 2 * sqrt_r * sqrt_lam_prod;
        bess_arg = bess_num ./ denom;

        abs_bess = abs(bess_arg);
        log_combined = exp_arg + abs_bess;
        bess_scaled = besseli(0, abs_bess, 1);

        prefactor = 2 * y_std * one_minus_rho2 ./ (1e4 * (1 + r).^2);
        f_std = prefactor .* exp(log_combined) .* bess_scaled;
        f(valid) = max(f_std / sig_ratio, 0);
    end

    f = reshape(f, sz);

end
