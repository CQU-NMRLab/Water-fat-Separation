function [seqParam] = read_seq_param(totalParam,paramName)
%READ_SEQ_PARAM �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    
iStart = strfind(totalParam,paramName);
if ~isempty(iStart)
    ik= iStart;
    while totalParam(ik)~=':'
        ik=ik+1;
    end
    tmp = totalParam(iStart:(ik-1));
    tmpParam = strsplit(tmp);
    seqParam = tmpParam{3};
else
    seqParam = ' ';
end
end

