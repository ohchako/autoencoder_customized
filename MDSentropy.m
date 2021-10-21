function [pvalue] = MDSentropy(R1_real_rate,Idx,x_test);


R1_real_rate_1 = R1_real_rate(:,1);
R1_real_rate_2 = R1_real_rate(:,2);
R1_real_rate_3 = R1_real_rate(:,3);
R1_real_rate_4 = R1_real_rate(:,4);
R1_real_rate_5 = R1_real_rate(:,5);

R1_real_rate_1_log = real(log2(R1_real_rate(:,1)));
R1_real_rate_2_log = real(log2(R1_real_rate(:,2)));
R1_real_rate_3_log = real(log2(R1_real_rate(:,3)));
R1_real_rate_4_log = real(log2(R1_real_rate(:,4)));
R1_real_rate_5_log = real(log2(R1_real_rate(:,5)));


% REAL entropy
clear ch_rate1 ch_rate1_sum ch_rate1_sum ch_rate1_sumsum P1 H1 H1_sum
ch_rate1 = R1_real_rate_1_log(Idx);
ch_rate1_sum = sum(ch_rate1,2);
ch_rate1_sumsum = sum(ch_rate1_sum);
P1 = ch_rate1_sum ./ ch_rate1_sumsum;
H1 = -(P1.*log2(P1));
H1 = real(H1);
H1_sum = sum(H1); % real entropy
% 
% clear ch_rate2 ch_rate2_sum ch_rate2_sum ch_rate2_sumsum P2 H2 H2_sum
% ch_rate2 = R1_real_rate_2_log(Idx);
% ch_rate2_sum = sum(ch_rate2,2);
% ch_rate2_sumsum = sum(ch_rate2_sum);
% P2 = ch_rate2_sum ./ ch_rate2_sumsum;
% H2 = -(P2.*log2(P2));
% H2 = real(H2);
% H2_sum = sum(H2); % real entropy

% clear ch_rate3 ch_rate3_sum ch_rate3_sum ch_rate3_sumsum P3 H3 H3_sum
% ch_rate3 = R1_real_rate_3_log(Idx);
% ch_rate3_sum = sum(ch_rate3,2);
% ch_rate3_sumsum = sum(ch_rate3_sum);
% P3 = ch_rate3_sum ./ ch_rate3_sumsum;
% H3 = -(P3.*log2(P3));
% H3 = real(H3);
% H3_sum = sum(H3); % real entropy
% 
% clear ch_rate4 ch_rate4_sum ch_rate4_sum ch_rate4_sumsum P4 H4 H4_sum
% ch_rate4 = R1_real_rate_4_log(Idx);
% ch_rate4_sum = sum(ch_rate4,2);
% ch_rate4_sumsum = sum(ch_rate4_sum);
% P4 = ch_rate4_sum ./ ch_rate4_sumsum;
% H4 = -(P4.*log2(P4));
% H4 = real(H4);
% H4_sum = sum(H4); % real entropy
% 
% clear ch_rate5 ch_rate5_sum ch_rate5_sum ch_rate5_sumsum P5 H5 H5_sum
% ch_rate5 = R1_real_rate_5_log(Idx);
% ch_rate5_sum = sum(ch_rate5,2);
% ch_rate5_sumsum = sum(ch_rate5_sum);
% P5 = ch_rate5_sum ./ ch_rate5_sumsum;
% H5 = -(P5.*log2(P5));
% H5 = real(H5);
% H5_sum = sum(H5); % real entropy


% SHUFFLE entropy
clear I
I = [1:size(x_test,1)]';
%IDX={};
shuffle_n = 1e6;

for j=1:shuffle_n
    for i=1:size(x_test,1)
        Idx_s(i,:) =randperm(size(x_test,1),4);
    end
    clear i
    Idx_s =[I Idx_s];

    for i=1:size(x_test,1)
    
        if Idx_s(i,1) == Idx_s(i,2)
            Idx_s(i,2) = randperm(size(x_test,1),1);
        elseif Idx_s(i,1) == Idx_s(i,3)
            Idx_s(i,3) = randperm(size(x_test,1),1);
        elseif Idx_s(i,1) == Idx_s(i,4)
            Idx_s(i,4) = randperm(size(x_test,1),1);
        elseif Idx_s(i,1) == Idx_s(i,5)
            Idx_s(i,5) = randperm(size(x_test,1),1);
        end
        
    end
%     idx = num2cell(Idx_s,1);
%     IDX=[IDX;idx];

ch_rate1_s = R1_real_rate_1_log(Idx_s);
ch_rate1_sum_s = sum(ch_rate1_s,2);
ch_rate1_sumsum_s = sum(ch_rate1_sum_s);
P1_s = ch_rate1_sum_s ./ ch_rate1_sumsum_s;
H1_s = -(P1_s.*log2(P1_s));
H1_s = real(H1_s);
H1_sum_s(j,:)=sum(H1_s); % shuffle entropy

% ch_rate2_s = R1_real_rate_2_log(Idx_s);
% ch_rate2_sum_s = sum(ch_rate2_s,2);
% ch_rate2_sumsum_s = sum(ch_rate2_sum_s);
% P2_s = ch_rate2_sum_s ./ ch_rate2_sumsum_s;
% H2_s = -(P2_s.*log2(P2_s));
% H2_s = real(H2_s);
% H2_sum_s(j,:)=sum(H2_s); % shuffle entropy
% 
% % ch_rate3_s = R1_real_rate_3_log(Idx_s);
% % ch_rate3_sum_s = sum(ch_rate3_s,2);
% % ch_rate3_sumsum_s = sum(ch_rate3_sum_s);
% % P3_s = ch_rate3_sum_s ./ ch_rate3_sumsum_s;
% % H3_s = -(P3_s.*log2(P3_s));
% % H3_s = real(H3_s);
% % H3_sum_s(j,:)=sum(H3_s); % shuffle entropy
% 
% ch_rate4_s = R1_real_rate_4_log(Idx_s);
% ch_rate4_sum_s = sum(ch_rate4_s,2);
% ch_rate4_sumsum_s = sum(ch_rate4_sum_s);
% P4_s = ch_rate4_sum_s ./ ch_rate4_sumsum_s;
% H4_s = -(P4_s.*log2(P4_s));
% H4_s = real(H4_s);
% H4_sum_s(j,:)=sum(H4_s); % shuffle entropy
% 
% ch_rate5_s = R1_real_rate_5_log(Idx_s);
% ch_rate5_sum_s = sum(ch_rate5_s,2);
% ch_rate5_sumsum_s = sum(ch_rate5_sum_s);
% P5_s = ch_rate5_sum_s ./ ch_rate5_sumsum_s;
% H5_s = -(P5_s.*log2(P5_s));
% H5_s = real(H5_s);
% H5_sum_s(j,:)=sum(H5_s); % shuffle entropy

clear ch_rate1_s ch_rate1_sum_s ch_rate1_sumsum_s P1_s H1_s ch_rate2_s ch_rate2_sum_s ch_rate2_sumsum_s P2_s H2_s ch_rate3_s ch_rate3_sum_s ch_rate3_sumsum_s P3_s H3_s
clear ch_rate4_s ch_rate4_sum_s ch_rate4_sumsum_s P4_s H4_s ch_rate5_s ch_rate5_sum_s ch_rate5_sumsum_s P5_s H5_s
clear Idx_s idx
end
clear j
clear entropy_s_sort1 entropy_s_sort2 entropy_s_sort3 entropy_s_sort4
entropy_s_sort1=sort(H1_sum_s);
% entropy_s_sort2=sort(H2_sum_s);
% entropy_s_sort3=sort(H3_sum_s);
% entropy_s_sort4=sort(H4_sum_s);
% entropy_s_sort5=sort(H5_sum_s);

clear entropy_s_inf95_1 entropy_s_inf95_2 entropy_s_inf95_3 entropy_s_inf95_4 entropy_s_sup95_1 entropy_s_sup95_2 entropy_s_sup95_3
entropy_s_inf95_1=entropy_s_sort1(shuffle_n*0.05+1);
% entropy_s_inf95_2=entropy_s_sort2(shuffle_n*0.05+1);
% % entropy_s_inf95_3=entropy_s_sort3(shuffle_n*0.05+1);
% entropy_s_inf95_4=entropy_s_sort4(shuffle_n*0.05+1);
% entropy_s_inf95_5=entropy_s_sort5(shuffle_n*0.05+1);

% ＜検定＞mmのリアルエントロピーをサロゲートの95%信頼区間の下限と比較する（）
% エントロピーの値が小さいほど、分散していない（95%の下限より小さくなっていてほしい！）
% if H1_sum < entropy_s_inf95_1
%     disp('SIGNIFICANT!!!');
% else
%     disp('Not significant...');
% end
% if H2_sum < entropy_s_inf95_2
%     disp('SIGNIFICANT!!!');
% else
%     disp('Not significant...');
% end
% if H3_sum < entropy_s_inf95_3
%     disp('SIGNIFICANT!!!');
% else
%     disp('Not significant...');
% end
% if H4_sum < entropy_s_inf95_4
%     disp('SIGNIFICANT!!!');
% else
%     disp('Not significant...');
% end
% if H5_sum < entropy_s_inf95_5
%     disp('SIGNIFICANT!!!');
% else
%     disp('Not significant...');
% end

clear surrogate_p1 surrogate_p2 surrogate_p3 surrogate_p4
surrogate_p1=entropy_s_sort1;
% surrogate_p2=entropy_s_sort2;
% % surrogate_p3=entropy_s_sort3;
% surrogate_p4=entropy_s_sort4;
% surrogate_p5=entropy_s_sort5;

clear Pv
Pv=(([1:shuffle_n]/shuffle_n)');
surrogate_p1(:,2)=Pv;
% surrogate_p2(:,2)=Pv;
% % surrogate_p3(:,2)=Pv;
% surrogate_p4(:,2)=Pv;
% surrogate_p5(:,2)=Pv;

clear IDX1 IDX2 IDX3 IDX4 IDX5
[IDX1]=knnsearch(surrogate_p1(:,1),H1_sum);
% [IDX2]=knnsearch(surrogate_p2(:,1),H2_sum);
% % [IDX3]=knnsearch(surrogate_p3(:,1),H3_sum);
% [IDX4]=knnsearch(surrogate_p4(:,1),H4_sum);
% [IDX5]=knnsearch(surrogate_p5(:,1),H5_sum);

clear ao1 ao2 ao3 ao4 ao5 fao1 fao2 fao3 fao4 fao5 mfao1 mfao2 mfao3 mfao4 mfao5 pvalue_all
ao1=round(surrogate_p1(IDX1),4);
% ao2=round(surrogate_p2(IDX2),4);
% % ao3=round(surrogate_p3(IDX3),4);
% ao4=round(surrogate_p4(IDX4),4);
% ao5=round(surrogate_p5(IDX5),4);

fao1=find(round(surrogate_p1(:,1),4)==ao1);
% fao2=find(round(surrogate_p2(:,1),4)==ao2);
% % fao3=find(round(surrogate_p3(:,1),4)==ao3);
% fao4=find(round(surrogate_p4(:,1),4)==ao4);
% fao5=find(round(surrogate_p5(:,1),4)==ao5);


mfao1=max(fao1);
% mfao2=max(fao2);
% mfao3=max(fao3);
% mfao4=max(fao4);
% mfao5=max(fao5);


pvalue = [Pv(mfao1)]';% Pv(mfao2)   Pv(mfao5)[
