function result = opendss_powerflow(dssFilePath, loadScaleFactors, DSSObj)
% OPENDSS_POWERFLOW 调用OpenDSS求解三相不平衡配电网潮流
%
% 语法:
%   result = opendss_powerflow(dssFilePath, loadScaleFactors, DSSObj)
%
% 输入:
%   dssFilePath      - OpenDSS主文件的完整路径 (string)
%   loadScaleFactors - 各负荷的有功/无功缩放因子 [nLoad x 2]
%                      第1列: 有功功率 kW
%                      第2列: 无功功率 kvar
%   DSSObj           - 已初始化的OpenDSS COM对象 (可选，传入以复用)
%
% 输出:
%   result - 结构体，包含:
%     .V_complex     - 各计算节点的复电压 [nCalcNode x 1] (复数, kV)
%     .V_pu          - 各计算节点的电压幅值 [nCalcNode x 1] (p.u.)
%     .busNames      - 母线名称 cell 数组
%     .nodeNames     - 计算节点名称 cell 数组 (busname.phase)
%     .nCalcNode     - 计算节点总数
%     .Y3ph          - 三相系统导纳矩阵 [nCalcNode x nCalcNode] (复数)
%     .nodeOrder     - OpenDSS系统Y矩阵的节点排列顺序
%     .lineNames     - 线路名称 cell 数组
%     .linePower     - 线路首末端功率 [nLine x 6]
%                      (P_from, Q_from, P_to, Q_to, P_loss, Q_loss) kW/kvar
%
% 版本: v1.0
% 日期: 2026-02-07

    arguments
        dssFilePath (1,1) string {mustBeNonempty}
        loadScaleFactors (:,2) double
        DSSObj = []
    end

    %% 1. 初始化 OpenDSS COM 接口
    if isempty(DSSObj)
        DSSObj = actxserver('OpenDSSEngine.DSS');
        assert(DSSObj.Start(0), '[OpenDSS] 启动失败');
    end

    DSSText    = DSSObj.Text;
    DSSCircuit = DSSObj.ActiveCircuit;

    %% 2. 编译电路文件
    DSSText.Command = ['Compile "', char(dssFilePath), '"'];

    DSSCircuit = DSSObj.ActiveCircuit;
    assert(~isempty(DSSCircuit), '[OpenDSS] 电路编译失败: %s', dssFilePath);

    %% 3. 设置负荷
    iLoad = DSSCircuit.Loads;
    nLoad = size(loadScaleFactors, 1);

    idx = iLoad.First;
    loadIdx = 1;
    while idx > 0 && loadIdx <= nLoad
        iLoad.kW   = loadScaleFactors(loadIdx, 1);
        iLoad.kvar = loadScaleFactors(loadIdx, 2);
        loadIdx = loadIdx + 1;
        idx = iLoad.Next;
    end

    %% 4. 求解潮流
    DSSSolution = DSSCircuit.Solution;
    DSSSolution.Solve;
    assert(DSSSolution.Converged, '[OpenDSS] 潮流不收敛');

    %% 5. 提取节点电压
    % AllBusVolts 返回所有节点的复电压 [Re1, Im1, Re2, Im2, ...]
    allVolts = DSSCircuit.AllBusVolts;
    % 重组为复数
    V_all = allVolts(1:2:end) + 1i * allVolts(2:2:end);

    % 获取节点名称和排列
    allNodeNames = DSSCircuit.AllNodeNames;
    nCalcNode = length(V_all);

    % 获取母线名称
    allBusNames = DSSCircuit.AllBusNames;

    result.V_complex  = V_all(:);
    result.V_pu       = DSSCircuit.AllBusVmagPu';
    result.busNames   = allBusNames;
    result.nodeNames  = allNodeNames;
    result.nCalcNode  = nCalcNode;

    %% 6. 提取系统导纳矩阵
    % SystemY 返回稀疏格式的Y矩阵 [Re1, Im1, Re2, Im2, ...]
    sysY = DSSCircuit.SystemY;
    Y_vec = sysY(1:2:end) + 1i * sysY(2:2:end);

    % 重组为方阵
    nY = round(sqrt(length(Y_vec)));
    Y3ph = reshape(Y_vec, [nY, nY]);

    % SystemY 的节点排列可能与 AllNodeNames 不同
    % 获取 Y 矩阵的节点顺序
    nodeOrderY = DSSCircuit.YNodeOrder;

    result.Y3ph      = Y3ph;
    result.nodeOrder  = nodeOrderY;

    %% 7. 提取线路功率
    iLine = DSSCircuit.Lines;
    nLine = iLine.Count;
    lineNames  = cell(nLine, 1);
    linePower  = zeros(nLine, 6); % Pf, Qf, Pt, Qt, Ploss, Qloss

    idx = iLine.First;
    lineIdx = 1;
    while idx > 0
        lineNames{lineIdx} = iLine.Name;
        % 需要通过 CktElement 获取功率
        DSSCircuit.SetActiveElement(['Line.' iLine.Name]);
        elem = DSSCircuit.ActiveCktElement;
        powers = elem.Powers; % [P1,Q1,P2,Q2,...] 各相功率

        % 总首端功率 = 前半部分之和, 总末端功率 = 后半部分之和
        nTermPhases = length(powers) / 2;
        halfIdx = nTermPhases / 2;

        P_from = sum(powers(1:2:nTermPhases));
        Q_from = sum(powers(2:2:nTermPhases));
        P_to   = sum(powers(nTermPhases+1:2:end));
        Q_to   = sum(powers(nTermPhases+2:2:end));

        linePower(lineIdx, :) = [P_from, Q_from, P_to, Q_to, ...
                                  P_from + P_to, Q_from + Q_to];
        lineIdx = lineIdx + 1;
        idx = iLine.Next;
    end

    result.lineNames = lineNames;
    result.linePower = linePower;
    result.nLine     = nLine;

end
