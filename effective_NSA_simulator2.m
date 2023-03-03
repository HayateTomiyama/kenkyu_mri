% n�摜�p����2���q�������s���ꍇ�̎���NSA�V�~�����[�V��������
function [effective_NSA]=effective_NSA_simulator2(sample,num)
% ���͎͂��̂悤�ɂȂ��Ă����
% sample�F0�`�΂����������Č��ʂ��o�͂��邩
% num�F�摜��

effective_NSA=zeros(1,sample);
% ����NSA�̏�����

% ����for���ł͎���NSA���Z�o���Ă����
for k=1:sample
    theta=k*pi/sample;
    effective_NSA(1,k)=num-(1-cos(num*theta))/(num*(1-cos(theta)));
end

end