% n画像用いて2分子分離を行う場合の実効NSAシミュレーションだよ
function [effective_NSA]=effective_NSA_simulator2(sample,num)
% 入力は次のようになっているよ
% sample：0～πを何等分して結果を出力するか
% num：画像数

effective_NSA=zeros(1,sample);
% 実効NSAの初期化

% 次のfor文では実効NSAを算出しているよ
for k=1:sample
    theta=k*pi/sample;
    effective_NSA(1,k)=num-(1-cos(num*theta))/(num*(1-cos(theta)));
end

end