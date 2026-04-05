function V_complex = eval_opendss_single(DSSObj, dssFile, x_vec, ...
    load_sequence, loadKW_base, loadKvar_base, nLoad)
% EVAL_OPENDSS_SINGLE 运行一次OpenDSS三相潮流，返回复电压向量
%
% 语法:
%   V_complex = eval_opendss_single(DSSObj, dssFile, x_vec, ...)
%
% 输入:
%   DSSObj        - 已初始化的OpenDSS COM对象
%   dssFile       - OpenDSS主文件路径 (char)
%   x_vec         - [1 x m] 或 [m x 1] 负荷缩放因子向量
%   load_sequence - [nLoad x 1] 各负荷对应的缩放因子列索引
%   loadKW_base   - [nLoad x 1] 各负荷基准有功功率 (kW)
%   loadKvar_base - [nLoad x 1] 各负荷基准无功功率 (kvar)
%   nLoad         - 负荷数量
%
% 输出:
%   V_complex - [nc x 1] 各计算节点复电压 (V)
%
% 版本: v1.0
% 日期: 2026-02-08

    arguments
        DSSObj
        dssFile       (1,:) char
        x_vec         (:,1) double
        load_sequence (:,1) double
        loadKW_base   (:,1) double
        loadKvar_base (:,1) double
        nLoad         (1,1) double {mustBePositive, mustBeInteger}
    end

    DSSText = DSSObj.Text;
    DSSText.Command = ['Compile "', dssFile, '"'];
    DSSCircuit = DSSObj.ActiveCircuit;

    iLoad = DSSCircuit.Loads;
    idx = iLoad.First;
    li = 1;
    while idx > 0 && li <= nLoad
        col = load_sequence(li);
        iLoad.kW   = loadKW_base(li) * x_vec(col);
        iLoad.kvar = loadKvar_base(li) * x_vec(col);
        li = li + 1;
        idx = iLoad.Next;
    end

    DSSCircuit.Solution.Solve;

    allVolts = DSSCircuit.AllBusVolts;
    V_complex = (allVolts(1:2:end) + 1i * allVolts(2:2:end)).';  % [nc x 1]

end
