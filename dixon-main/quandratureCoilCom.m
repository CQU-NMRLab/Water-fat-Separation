function [comData] = quandratureCoilCom(coil1data,coil2data,absCeo)
%% QuandratureCoilCom 正交线圈数据组合叠加
%   输入参数：
%       coil1data       线圈1的k空间数据
%       coil2data       线圈2的k空间数据
%       absCeo          叠加系数幅值，角度采用quandratureCoilCeo函数结果
%   输出参数：
%       comCeo          线圈1/2之间的拟合系数，comCeo*coil1=coil2

    if nargin==2
        absCeo = 1;
    end
    coe = quandratureCoilCoe(coil1data,coil2data);  % 计算拟合系数
    newCoe = absCeo*exp(sqrt(-1)*angle(coe));
    comData = coil1data+newCoe.*coil2data;
end

