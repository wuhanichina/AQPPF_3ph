function f = gx2pdf_davies(x,w,k,lambda,m,s,varargin)  
% GX2PDF_DAVIES 计算广义卡方分布的 PDF（使用 Davies 方法）
% 这是一个加权非中心卡方分布和正态分布之和的 PDF
%
% Required inputs:
%   x         评估点，或 'full' 返回完整范围的 pdf
%   w         非中心卡方的权重行向量
%   k         非中心卡方的自由度行向量
%   lambda    非中心卡方的非中心参数行向量（均值平方和）
%   m         正态项的均值
%   s         正态项的标准差
%
% Optional name-value inputs:
%   method    'diff' (default) 使用 CDF 数值微分
%             'conv' 使用非中心卡方 PDF 卷积
%   dx        步长精度（用于卷积或微分）
%   AbsTol    使用 'diff' 方法时的绝对误差容限
%   RelTol    使用 'diff' 方法时的相对误差容限
%
% Output:
%   f         计算得到的 PDF 值
%
% See also:
%   gx2cdf, gx2cdf_davies, gx2cdf_imhof, gx2cdf_ruben

    parser = inputParser;
    addRequired(parser,'x',@(x) isreal(x) || strcmpi(x,'full'));
    addRequired(parser,'w',@(x) isreal(x) && isrow(x));
    addRequired(parser,'k',@(x) isreal(x) && isrow(x));
    addRequired(parser,'lambda',@(x) isreal(x) && isrow(x));
    addRequired(parser,'m',@(x) isreal(x) && isscalar(x));
    addRequired(parser,'s',@(x) isreal(x) && isscalar(x));
    addParameter(parser,'method','diff',@(s) strcmpi(s,'diff') || strcmpi(s,'conv'));
    addParameter(parser,'AbsTol',1e-10,@(x) isreal(x) && isscalar(x) && (x>=0));
    addParameter(parser,'RelTol',1e-6,@(x) isreal(x) && isscalar(x) && (x>=0));
    
    % 获取分布的方差以确定默认步长
    [~,v] = gx2stat(w,k,lambda,m,s);
    addParameter(parser,'dx',sqrt(v)/1e2,@(x) isreal(x) && isscalar(x) && (x>=0)); % 默认步长为 std/100
    
    parse(parser,x,w,k,lambda,m,s,varargin{:});
    dx = parser.Results.dx;
    
    % 使用数值微分计算 PDF
    p_left = gx2cdf_davies(x-dx,w,k,lambda,m,s);
    p_right = gx2cdf_davies(x+dx,w,k,lambda,m,s);
    f = (p_right-p_left)/(2*dx);
    f = max(f,0);
end

