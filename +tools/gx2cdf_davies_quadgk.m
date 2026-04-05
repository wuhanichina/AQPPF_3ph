function p = gx2cdf_davies_quadgk(x, w, k, lambda, m, s, opts)
% GX2CDF_DAVIES_QUADGK 计算广义卡方分布CDF（Davies法，quadgk积分版）
%
% 目的：避免 MATLAB integral() 在某些版本/参数下触发“区间数上限”警告，
%      同时不修改 gncs/gx2cdf_davies.m。
%
% 模型：X = sum_i w_i * ChiSq(k_i, lambda_i) + N(m, s^2)
% 返回：p = Pr(X <= x)
%
% 输入:
%   x      - 标量或向量（实数）
%   w,k,lambda - 行向量参数（长度一致）
%   m,s    - 标量
%   opts   - 可选结构体（第7个位置参数）
%            .side: 'lower' 或 'upper'（默认'lower'）
%            .AbsTol: quadgk绝对误差（默认1e-10）
%            .RelTol: quadgk相对误差（默认1e-6）
%            .MaxIntervalCount: 最大区间数（默认200000）
%
% 输出:
%   p      - 与x同形状
%
% 版本: v1.0
% 日期: 2025-12-20

arguments
    x (:,1) double {mustBeReal}
    w (1,:) double {mustBeReal}
    k (1,:) double {mustBeReal}
    lambda (1,:) double {mustBeReal}
    m (1,1) double {mustBeReal}
    s (1,1) double {mustBeReal, mustBeNonnegative}
    opts (1,1) struct = struct()
end

% 默认选项（Fail-Fast：缺字段即用默认，不做 try-catch）
if ~isfield(opts, 'side'); opts.side = 'lower'; end
if ~isfield(opts, 'AbsTol'); opts.AbsTol = 1e-10; end
if ~isfield(opts, 'RelTol'); opts.RelTol = 1e-6; end
if ~isfield(opts, 'MaxIntervalCount'); opts.MaxIntervalCount = 200000; end

assert(ischar(opts.side) || isstring(opts.side), '[gx2cdf_davies_quadgk] opts.side 类型错误');
opts.side = char(lower(string(opts.side)));
assert(any(strcmp(opts.side, {'lower','upper'})), '[gx2cdf_davies_quadgk] opts.side 必须为 lower/upper');
assert(isscalar(opts.AbsTol) && opts.AbsTol >= 0, '[gx2cdf_davies_quadgk] opts.AbsTol 必须为非负标量');
assert(isscalar(opts.RelTol) && opts.RelTol >= 0, '[gx2cdf_davies_quadgk] opts.RelTol 必须为非负标量');
assert(isscalar(opts.MaxIntervalCount) && opts.MaxIntervalCount > 0, '[gx2cdf_davies_quadgk] opts.MaxIntervalCount 必须为正标量');

assert(isrow(w) && isrow(k) && isrow(lambda), '[gx2cdf_davies_quadgk] w/k/lambda 必须为行向量');
assert(numel(w) == numel(k) && numel(w) == numel(lambda), '[gx2cdf_davies_quadgk] 参数长度不一致');

w = w(:); k = k(:); lambda = lambda(:); % 列向量

% Davies integrand
    function f = davies_integrand(u, x0)
        % 避免 u=0 处 1/u 的数值问题：quadgk 可能会采样到 u=0
        u_safe = max(u, realmin);
        theta = sum(k .* atan(w .* u_safe) + (lambda .* (w .* u_safe)) ./ (1 + (w.^2) .* (u_safe.^2)), 1) / 2 - u_safe .* x0 / 2;
        rho = prod(((1 + (w.^2) .* (u_safe.^2)).^(k/4)) .* exp(((w.^2) .* (u_safe.^2) .* lambda) ./ (2 * (1 + (w.^2) .* (u_safe.^2)))), 1) .* exp((u_safe.^2) * s^2 / 8);
        f = sin(theta) ./ (u_safe .* rho);
    end

% 逐点积分
p = zeros(size(x));
for ii = 1:numel(x)
    x0 = x(ii) - m;
    I = quadgk(@(u) davies_integrand(u, x0), 0, inf, ...
        'AbsTol', opts.AbsTol, 'RelTol', opts.RelTol, 'MaxIntervalCount', opts.MaxIntervalCount);
    if strcmpi(opts.side, 'lower')
        p(ii) = 0.5 - I / pi;
    else
        p(ii) = 0.5 + I / pi;
    end
end

% 夹到[0,1]并fail-fast检查
assert(all(isfinite(p)), '[gx2cdf_davies_quadgk] 输出包含非有限值');
p = max(p, 0);
p = min(p, 1);

end



