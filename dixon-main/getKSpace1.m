function [Image_pad,seqParam] = getKSpace1(mrdFile)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
fid = fopen(mrdFile,'r');               %定义文件ID
%读文件头信息（126 byte）
Dimenstion = fread(fid,4,'int32');        %Dimenstion为MRD文件前四个数据
no_samples = Dimenstion(1);               %number of data samples
no_views =  Dimenstion(2);                %number of views (2D data)
no_views_2= Dimenstion(3);                %number of secondary views (3D data)
no_slices = Dimenstion(4);                %number of slices
fseek(fid,18,'bof');                      %从文件开始(Beginning of file)偏移18个字节的位置
datatype=fread(fid,1, 'uint16');          %数据类型代码
datatype = dec2hex(datatype);             %转换为16进制
fseek(fid,152,'bof');                     %从文件开始偏移152个字节的位置
Dimenstion(5:6) = fread(fid,2, 'int32');  %MRD文件Dimenstion的五和六的数据
no_echoes =  Dimenstion(5);               %number of echoes
no_experiments =  Dimenstion(6);          %number of experiments
%读文件描述（126 byte）
fseek(fid,256,'bof');                     %从文件开始偏移256个字节的位置
text=fread(fid,256);                      %256 bytes ASCII Text, 0 terminated string
%确定数据类型代码
if size(datatype,2)>1                     
    datatypecde = datatype(2);             
    iscomplex = 2;                         %是否为复数：2为复数
else
    datatypecde = datatype(1); 
    iscomplex = 1;
end
%确定数据的类型
switch datatypecde
    case '0'
        dataformat = 'uchar';   datasize = 1; % size in bytes
    case '1'
        dataformat = 'schar';   datasize = 1; % size in bytes
    case '2'
        dataformat = 'short';   datasize = 2; % size in bytes
    case '3'
        dataformat = 'int16';   datasize = 2; % size in bytes √
    case '4'
        dataformat = 'int32';   datasize = 4; % size in bytes
    case '5'
        dataformat = 'float32'; datasize = 4; % size in bytes
    case '6'
        dataformat = 'double';  datasize = 8; % size in bytes
end
%计算数据的数量
no_data = no_experiments*no_echoes*no_slices*no_views_2*no_views*no_samples*iscomplex; %*datasize ;
[data_total, count] = fread(fid,no_data,dataformat);      %读入K空间数据
%% 获取所有参数内容
totalParam = fread(fid,[1,inf]);
totalParam(totalParam==10) = [];
totalParam(totalParam==13) = [];
totalParam = char(totalParam);
iStart = strfind(char(totalParam),'FOV ');
ik= iStart;
while totalParam(ik)~=':'
    ik=ik+1;
end
tmp = totalParam(iStart:(ik-1));
tmpParam = strsplit(tmp);
seqParam.FOV = tmpParam{2};

seqParam.te = read_seq_param(totalParam,'VAR te,');
seqParam.tr = read_seq_param(totalParam,'VAR tr,');
seqParam.tRgs = read_seq_param(totalParam,'VAR tRgs,');
seqParam.FA = read_seq_param(totalParam,'VAR alpha,');
seqParam.res = [num2str(no_samples),'*',num2str(no_views)];

fclose(fid);   %关闭文件

%判断数据点是否对应
if (count~=no_data)
    h = msgbox('The data have a problem...'); 
end
%对数据进行整理
if iscomplex == 2
    a=1:no_data/2;            
    data_real = data_total(2*a-1);    %data_real(1.3.5.7.9...)
    
    data_imag = data_total(2*a);      %data_imag(2.4.6.8...)
%     data_imag = medfilt1(data_imag,5);
%     data_real = medfilt1(data_real,5);
    data_C = data_real+data_imag*i;
else
    data_C = data_total;
end
%构成K空间数据的矩阵
n = 0;
for k = 1:no_experiments
    for j = 1:no_echoes
        for q=1:no_slices
            for a = 1:no_views
                for b = 1:no_views_2
                    for c = 1:no_samples
                        n = n + 1;
                        Kspacedata(a,c,b,q,j,k) = data_C(n);
                    end
                end
            end
        end
    end
end
t=1;
for k = 1:no_experiments
    for j = 1:no_echoes
        for q =1:no_slices
            for m = 1: no_views_2
                KspaceData(:,:,t)=Kspacedata(:,:,m,q,j,k);
                t=t+1;
            end
        end
    end
end
DIM1=512;
DIM2=512;
DIM1_pad=round(DIM1/2-no_views/2);
DIM2_pad=round(DIM2/2-no_samples/2);
DIM3_pad=0;
Image_pad=zeros(DIM1,DIM2,no_views_2);
for w = 1: no_views_2*no_slices*no_echoes*no_experiments
    Image_pad(DIM1_pad:DIM1_pad+no_views-1,DIM2_pad:DIM2_pad+no_samples-1,w)=KspaceData(:,:,w);
end
end

