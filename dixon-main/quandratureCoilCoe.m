function [comCeo] = quandratureCoilCoe(coil1data,coil2data)
%% quandratureCoilCeo 根据磁共振信号求解两个线圈数据之间的拟合系数
%       采用磁共振信号数据中磁共振信号最大的数据部分进行拟合。
%   输入参数：
%       coil1data       线圈1的k空间数据
%       coil2data       线圈2的k空间数据
%   输出参数：
%       comCeo          线圈1/2之间的拟合系数，comCeo*coil1=coil2
    tmp = size(coil1data,3);
    if tmp == 1
        peInx = floor(size(coil1data,1)/2)+1;
        sig1 = coil1data(peInx,:);
        sig2 = coil2data(peInx,:);
        comCeo = sig1'\sig2';
    else
        peInx = floor(size(coil1data,1)/2)+1;
        seInx = floor(size(coil1data,3)/2)+1;
        sig1 = coil1data(peInx,:,seInx);
        sig2 = coil2data(peInx,:,seInx);
        comCeo = sig1'\sig2';
    end

end

