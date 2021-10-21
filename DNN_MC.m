%% DNN_MC.m
%% usage; peak50ms: N x 2001 matrix(for vitro),N x 4001 matrix(for vivo)  
%% update 211021
%% modulate LFP trace → as x_test(after processing python)
% scaling in the range of [0-1]（210311）

peak50ms_n=normalize(peak50ms,2,'range');

csvwrite('slice2_5cells_102.txt',peak50ms_n);% <---- modify


%% modulate vm trace

peak50ms_vm1=peak50ms_vm(:,:,1);
peak50ms_vm2=peak50ms_vm(:,:,2);
peak50ms_vm3=peak50ms_vm(:,:,3);
peak50ms_vm4=peak50ms_vm(:,:,4);
peak50ms_vm5=peak50ms_vm(:,:,5);

peak50ms_vm1_n=normalize(peak50ms_vm1,2,'range');
peak50ms_vm2_n=normalize(peak50ms_vm2,2,'range');
peak50ms_vm3_n=normalize(peak50ms_vm3,2,'range');
peak50ms_vm4_n=normalize(peak50ms_vm4,2,'range');
peak50ms_vm5_n=normalize(peak50ms_vm5,2,'range');

peak50ms_vm_n=cat(3,peak50ms_vm1_n,peak50ms_vm2_n,peak50ms_vm3_n,peak50ms_vm4_n,peak50ms_vm5_n);
csvwrite('slice2_vmtrace_normalized_102.txt',peak50ms_vm_n); % <---- modify


%% create random numbers for shuffle
for j=1:size(peak50ms_vm,3)
    for i=1:100
        R(:,i,j)=randperm(size(peak50ms_vm,1));
    end
end
csvwrite('slice2_Vmshuffleno_102.txt',R); % <---- modify

%% load output waveform learning with python→x_test

csvfiles=dir('*.csv');% be careful for the PATH!!
numfiles=length(csvfiles);
X_test=cell(1,numfiles);

for k=1:numfiles
    X_test{k}=importdata(csvfiles(k).name);
end
clear k numfiles

%% processing after training with python ## k-fold cross validation ##
% 210630
% Decoded_vm; output SW waveform
% X_test; original SW waveform (scaled 0-1)
% above two variables stored as 1×k cell（k = 5 in this cae）

% Average of the bottom 10-20% of LFP source data（①）
for j = 1 : size(X_test,2)
    for i = 1 : size(X_test{1,1},1)
        LFP_sort{:,j}(i,:) = sort(X_test{:,j}(i,:));
        LFP_sort_10_20{:,j}(i,:) = LFP_sort{:,j}(i,200:400);
        LFP_sort_mean{:,j}(i,:) = mean(LFP_sort_10_20{:,j}(i,:),2);
    end
end
clear i j

% Average of the bottom 10-20% of decoded SW waveform（②）
for j = 1 : size(Decoded_vm,2)
    for i = 1 : size(Decoded_vm{1,1},1)
        vmsort{:,j}(i,:) = sort(Decoded_vm{:,j}(i,:));
        vmsort_10_20{:,j}(i,:) = vmsort{:,j}(i,200:400);
        vmsort_mean{:,j}(i,:)= mean(vmsort_10_20{:,j}(i,:),2);
    end
end
clear i j

% correct baseline（② - ①）
vmsort_double = cell2mat(vmsort_mean);
LFPsort_double = cell2mat(LFP_sort_mean);
base_corr_comb = vmsort_double - LFPsort_double;

% substract the corrected baseline from the decoded SW waveform
for i = 1 : size(Decoded_vm,2)
    m_LFP(:,i) = {Decoded_vm{1,i} - base_corr_comb(:,i)};
end
clear i

% RMSE（original SW trace used as training data vs corrected decoded SW waveform）
for j=1:size(Decoded_vm,2)
    for i=1:size(X_test{1,1},1)
        RMSE_corr_cv{:,j}(i,:) = sqrt(immse(X_test{:,j}(i,:),m_LFP{:,j}(i,:)));
    end
end
clear i j
RMSE_corr_cv = cell2mat(RMSE_corr_cv);

RMSE_corr_cv_m = mean(RMSE_corr_cv,2);

%% For each combination、after shuffling ## import data ##
% calculate RMSE in each combination with correcting baseling
% usage; decoded SW waveform: size(x_test,1)×2001×N（or 1×N cell）
% 上述のLFPsort_mean_2は固定される
clear csvfiles numfiles

% import data
csvfiles=dir('*.csv');%% be careful for the PATH!!
numfiles=length(csvfiles);
decoded_vm=cell(1,numfiles);

for k=1:numfiles
    decoded_vm{k}=importdata(csvfiles(k).name);
end
clear k numfiles
%% For each combination、after shuffling　## import data ##
% calculate RMSE in each combination with correcting baseling
% usage; decoded SW waveform: size(x_test,1)×2001×N（or 1×N cell）
% 上述のLFPsort_mean_2は固定される
clear csvfiles numfiles

% import data
csvfiles=dir('*.csv');% be careful for the PATH!!
numfiles=length(csvfiles);
x_test=cell(1,numfiles);

for k=1:numfiles
    x_test{k}=importdata(csvfiles(k).name);
end
clear k numfiles
%% For each combination、after shuffling　## correct baseline & calculate RMSE ##
% ## k-fold cross-validation ##
% 5-fold cv、10-times shuffle each in python（total 50 of x_test and decoded_vm_s）
clear vmsort vmsort_10_20 vmsort_mean csvfiles numfiles vmsort_double base_corr_comb m_LFP_s LFPsort_double
clear LFP_sort LFP_sort_10_20 LFP_sort_mean combset shufflenum

combset = 4; % number of combinations
shufflenum = 10; % 1 cell prediction→100、other→10

% Average of the bottom 10-20% of LFP source data（①）
for j = 1 : size(x_test,2)
    for i = 1 : size(x_test{1,1},1)
        LFP_sort{:,j}(i,:) = sort(x_test{:,j}(i,:));
        LFP_sort_10_20{:,j}(i,:) = LFP_sort{:,j}(i,200:400);
        LFP_sort_mean{:,j}(i,:) = mean(LFP_sort_10_20{:,j}(i,:),2);
    end
end
clear i j

% Average of the bottom 10-20% of decoded SW waveform（②）
for j = 1 : size(decoded_vm,2)
    for i = 1 : size(decoded_vm{1,1},1)
        vmsort{:,j}(i,:) = sort(decoded_vm{:,j}(i,:));
        vmsort_10_20{:,j}(i,:) = vmsort{:,j}(i,200:400);
        vmsort_mean{:,j}(i,:)= mean(vmsort_10_20{:,j}(i,:),2);
    end
end
clear i j

% correct baseline（② - ①）
vmsort_double = cell2mat(vmsort_mean);
LFPsort_double = cell2mat(LFP_sort_mean);
base_corr_comb = vmsort_double - LFPsort_double;

% substract the corrected baseline from the decoded SW waveform
for i = 1 : size(decoded_vm,2)
    m_LFP_s(:,i) = {decoded_vm{1,i} - base_corr_comb(:,i)};
end
clear i

% RMSE（original SW trace used as training data vs corrected decoded SW waveform）
for j=1:size(decoded_vm,2)
    for i=1:size(x_test{1,1},1)
        RMSE_corr_cv_s{:,j}(i,:) = sqrt(immse(x_test{:,j}(i,:),m_LFP_s{:,j}(i,:)));
    end
end
clear i j

RMSE_corr_cv_s = cell2mat(RMSE_corr_cv_s);
RMSE_corr_cv_s = reshape(RMSE_corr_cv_s,[size(base_corr_comb,1)*5,combset*shufflenum]);

RMSE_corr_cv_s_m = mean(RMSE_corr_cv_s,2);

%% For each combination、after shuffling　## correct baseline & calculate RMSE ##
% #### shuffle-cvは取り込み方に注意 ####
% ## k-fold cross-validation ##
%  5-fold cv、10-times shuffle each in python（total 50 of x_test and decoded_vm_s）
clear vmsort vmsort_10_20 vmsort_mean csvfiles numfiles vmsort_double base_corr_comb m_LFP_s LFPsort_double
clear LFP_sort LFP_sort_10_20 LFP_sort_mean combset shufflenum
% 
% combset = 1; % number of combinations
% shufflenum = 10; % 1 cell predictionのときは100、それ以外は10！

% Average of the bottom 10-20% of LFP source data（①）
for j = 1 : size(x_test,2)
    for i = 1 : size(x_test{1,1},1)
        LFP_sort{:,j}(i,:) = sort(x_test{:,j}(i,:));
        LFP_sort_10_20{:,j}(i,:) = LFP_sort{:,j}(i,200:400);
        LFP_sort_mean{:,j}(i,:) = mean(LFP_sort_10_20{:,j}(i,:),2);
    end
end
clear i j

% Average of the bottom 10-20% of decoded SW waveform（②）
for j = 1 : size(decoded_vm,2)
    for i = 1 : size(decoded_vm{1,1},1)
        vmsort{:,j}(i,:) = sort(decoded_vm{:,j}(i,:));
        vmsort_10_20{:,j}(i,:) = vmsort{:,j}(i,200:400);
        vmsort_mean{:,j}(i,:)= mean(vmsort_10_20{:,j}(i,:),2);
    end
end
clear i j

% correct baseline（② - ①）
vmsort_double = cell2mat(vmsort_mean);
LFPsort_double = cell2mat(LFP_sort_mean);
base_corr_comb = vmsort_double - LFPsort_double;

% substract the corrected baseline from the decoded SW waveform
for i = 1 : size(decoded_vm,2)
    m_LFP_s(:,i) = {decoded_vm{1,i} - base_corr_comb(:,i)};
end
clear i

% RMSE（original SW trace used as training data vs corrected decoded SW waveform）
for j=1:size(decoded_vm,2)
    for i=1:size(x_test{1,1},1)
        RMSE_corr_cv_s{:,j}(i,:) = sqrt(immse(x_test{:,j}(i,:),m_LFP_s{:,j}(i,:)));
    end
end
clear i j

RMSE_corr_cv_s = cell2mat(RMSE_corr_cv_s);

for i = 1:10
    % RMSE_corr_cv_Re(:,1)が各cvにおける（つまり全SW）シャッフル1回目
    RMSE_corr_cv_Re(:,i) = [RMSE_corr_cv_s(:,i);RMSE_corr_cv_s(:,i+10);RMSE_corr_cv_s(:,i+20);RMSE_corr_cv_s(:,i+30);RMSE_corr_cv_s(:,i+40)];
end
clear i
clear RMSE_corr_cv_s
RMSE_corr_cv_s = RMSE_corr_cv_Re;
clear RMSE_corr_cv_Re

RMSE_corr_cv_s_m = mean(RMSE_corr_cv_s,2);

%% For each combination、after shuffling　## correct baseline & calculate RMSE ##
% ## k-fold cross-validation ##
% 5-fold cv、10-times shuffle each in python（total 50 of x_test and decoded_vm_s）
clear vmsort vmsort_10_20 vmsort_mean csvfiles numfiles vmsort_double base_corr_comb m_LFP_s LFPsort_double
clear LFP_sort LFP_sort_10_20 LFP_sort_mean

combset = 1;

% Average of the bottom 10-20% of LFP source data（①）
for j = 1 : size(x_test,2)
    for i = 1 : size(x_test{1,1},1)
        LFP_sort{:,j}(i,:) = sort(x_test{:,j}(i,:));
        LFP_sort_10_20{:,j}(i,:) = LFP_sort{:,j}(i,200:400);
        LFP_sort_mean{:,j}(i,:) = mean(LFP_sort_10_20{:,j}(i,:),2);
    end
end
clear i j

% Average of the bottom 10-20% of decoded SW waveform（②）
for j = 1 : size(decoded_vm,2)
    for i = 1 : size(decoded_vm{1,1},1)
        vmsort{:,j}(i,:) = sort(decoded_vm{:,j}(i,:));
        vmsort_10_20{:,j}(i,:) = vmsort{:,j}(i,200:400);
        vmsort_mean{:,j}(i,:)= mean(vmsort_10_20{:,j}(i,:),2);
    end
end
clear i j

% correct baseline（② - ①）
vmsort_double = cell2mat(vmsort_mean);
LFPsort_double = cell2mat(LFP_sort_mean);
base_corr_comb = vmsort_double - LFPsort_double;

% substract the corrected baseline from the decoded SW waveform
for i = 1 : size(decoded_vm,2)
    m_LFP_s(:,i) = {decoded_vm{1,i} - base_corr_comb(:,i)};
end
clear i


% RMSE（original SW trace used as training data vs corrected decoded SW waveform）
for j=1:size(decoded_vm,2)
    for i=1:size(x_test{1,1},1)
        RMSE_corr_cv{:,j}(i,:) = sqrt(immse(x_test{:,j}(i,:),m_LFP_s{:,j}(i,:)));
    end
end
clear i j


RMSE_corr_cv = cell2mat(RMSE_corr_cv);
RMSE_corr_cv = reshape(RMSE_corr_cv,[size(base_corr_comb,1)*5,combset]);

% RMSE_corr_cv_m = mean(RMSE_corr_cv,2);


%% cumulative probability (fig.3c)

% All = [RMSE_corr_cv];
All = [All;RMSE_corr_cv];
RMSE_corr_cv_s_re = reshape(RMSE_corr_cv_s,[size(RMSE_corr_cv_s,1)*size(RMSE_corr_cv_s,2),1]);
% All_s = [RMSE_corr_cv_s_re];
All_s = [All_s;RMSE_corr_cv_s_re];

clearvars -except All All_s


clear f x F X h p k
[f,x] = ecdf(All);
[F,X] = ecdf(All_s);
[h,p,k] = kstest2(x,X)
% 
figure;plot(x,f);hold on;plot(X,F);
ylabel('cumulative probability');xlabel('RMSE');title('all slices');xlim([0 0.25]);
legend('real','shuffle','Location','southeast');



%% fig2d (211020;fig2d→fig3d)
% 210715
% load '210715_fig2d.mat'

sz=100;
yyaxis right
semilogy(RMSE(:,1),RMSE(:,3),'color',[255/255 70/255 50/255]);hold on;scatter(RMSE(:,1),RMSE(:,3),sz,'.','MarkerEdgeColor',[255/255 70/255 50/255],'MarkerFaceColor',[255/255 70/255 50/255]);hold on;plot(RMSE(:,5),RMSE(:,4),'--','color',[255/255,70/255,50/255]);set(gca,'YDir','reverse');xlim([0.5 5.5]);ylim([10^(-50) 10^5]);
yyaxis left
plot(RMSE(:,1),RMSE(:,2),'k');hold on;scatter(RMSE(:,1),RMSE(:,2),sz,'k','.');ylim([0 0.15])


%% for 3 cells analysis
% update; 210331

% % 1細胞のときのRMSE(R1)と寄与率の比較（相関関係をみたい）
% % 1細胞のRMSEをチャンスレベル（シャッフル）を考慮して計算する
% % 各細胞のRMSEシャッフルの平均を計算する
% 3cells
shuffle_mean = [mean(RMSE_corr_3c1_s(:,1:100),2) mean(RMSE_corr_3c1_s(:,101:200),2) mean(RMSE_corr_3c1_s(:,201:300),2)];

% 4cells
shuffle_mean = [mean(RMSE_corr_4c1_s(:,1:100),2) mean(RMSE_corr_4c1_s(:,101:200),2) mean(RMSE_corr_4c1_s(:,201:300),2) mean(RMSE_corr_4c1_s(:,301:400),2)];

% 5cells
shuffle_mean = [mean(RMSE_corr_5c1_s(:,1:100),2) mean(RMSE_corr_5c1_s(:,101:200),2) mean(RMSE_corr_5c1_s(:,201:300),2) mean(RMSE_corr_5c1_s(:,301:400),2) mean(RMSE_corr_5c1_s(:,401:500),2)];


%★(S1-R1)/S1
for i = 1:size(shuffle_mean,1)
    for j= 1:size(shuffle_mean,2)
        R1_real(i,j) = shuffle_mean(i,j) - RMSE_corr_5c1(i,j); % <---- modify
        R1_real_rate(i,j) = R1_real(i,j)/shuffle_mean(i,j);
    end
end

% --- cosine similarity ---(fig.4c)
v=1:1:5; C=nchoosek(v,2);
for i=1:size(C,1)
    cs(i,:) = getCosineSimilarity(R1_real_rate(:,C(i,1)),R1_real_rate(:,C(i,2)));
end

cs_all = [cs_all;cs];

clearvars -except cs_all

figure;
histogram(cs_all,8);xlabel('Cosine Similarity');ylabel('Count');xlim([-1 1]);

cs_all_cosd = acosd(cs_all); % cosine similarity to angle(deg)
cs_all_cosd_r = deg2rad(cs_all_cosd);% angle(deg) to radian

[pval,v] = circ_vtest(cs_all_cosd_r,circ_ang2rad(90)) % modulated calculating method for pval

figure;
polarhistogram(cs_all_cosd_r,15,'FaceColor',[20/255,180/255,20/255]);thetalim([0 180]);


%% calculating D value 
% modify in each combinations
%210615
% update 210804

RMSE_corr_cv_3c1_s1=RMSE_corr_cv_3c1_s(:,1:100);
RMSE_corr_cv_3c1_s2=RMSE_corr_cv_3c1_s(:,101:200);
RMSE_corr_cv_3c1_s3=RMSE_corr_cv_3c1_s(:,201:300);

RMSE_corr_cv_3c1_s1_sort=sort(RMSE_corr_cv_3c1_s1,2);
RMSE_corr_cv_3c1_s2_sort=sort(RMSE_corr_cv_3c1_s2,2);
RMSE_corr_cv_3c1_s3_sort=sort(RMSE_corr_cv_3c1_s3,2);

RMSE_corr_3c1_s_sort={RMSE_corr_cv_3c1_s1_sort RMSE_corr_cv_3c1_s2_sort RMSE_corr_cv_3c1_s3_sort};

% KStest
RMSE_shuffle = {reshape(RMSE_corr_cv_3c1_s1,[],1) reshape(RMSE_corr_cv_3c1_s2,[],1) reshape(RMSE_corr_cv_3c1_s3,[],1)};

clear f x F X h p k
for i=1:size(RMSE_corr_cv_3c1,2)
    [f(:,i),x(:,i)] = ecdf(RMSE_corr_cv_3c1(:,i));
    [F(:,i),X(:,i)] = ecdf(RMSE_shuffle{1,i});
    [~,p(i,:),k(i,:)] = kstest2(x(:,i),X(:,i)); % pval,dval
end

%% hilus map (figS4A)

% update;210818
sizex=[1:1:87]';
h=imagesc('CData',sizex);
CData=h.CData;
cmap = jet(256);
cmin=min(CData(:));
cmax=max(CData(:));
m=length(cmap);
index = fix((CData-cmin)/(cmax-cmin)*m)+1;
RGB=ind2rgb(index,cmap);
RGB_stack = [RGB(:,:,1) RGB(:,:,2) RGB(:,:,3)];
% RGB_stack_c =RGB_stack .* 255;% for PP

% for check
% x=1:1:87;
% y=1:1:87;
% figure;gscatter(x,y,x,RGB_stack);

% import D value
[~,ind] = sort(D);
R = deg2rad(norm_angle(ind));
rho = norm_dist(ind);

% figS4A
figure;
for i=1:size(R,1)
pl = polarplot(R(i,:),rho(i,:),'.');thetalim([0 180]);rlim([0 1]);hold on;
pl.Color = RGB_stack(i,:);
pl.MarkerSize = 12;
end

RR=[deg2rad(norm_angle) norm_dist];

% figS4B --spatial entropy
[Idx,~] = knnsearch(RR,RR,'K',5);% get 4 points around, 5 points including own point.default=euclidean distance

D_val = D(Idx);
ch_rate1_sum = sum(D_val,2);
D_val_sumsum = sum(ch_rate1_sum);
P = ch_rate1_sum ./ D_val_sumsum;
H = -(P.*log2(P));
H_sum = sum(H); % real entropy


% shuffle Idx（N×5）（N is number of x_test）the most left side is fixed. Strictly shuffle of N x 4
I = [1:size(D,1)]';
for j=1:1e4
    for i=1:size(D,1)
        Idx_s(i,:) =randperm(size(D,1),4);
    end
    clear i
    Idx_s =[I Idx_s];

    for i=1:size(D,1)
    
        if Idx_s(i,1) == Idx_s(i,2)
            Idx_s(i,2) = randperm(size(D,1),1);
        elseif Idx_s(i,1) == Idx_s(i,3)
            Idx_s(i,3) = randperm(size(D,1),1);
        elseif Idx_s(i,1) == Idx_s(i,4)
            Idx_s(i,4) = randperm(size(D,1),1);
        elseif Idx_s(i,1) == Idx_s(i,5)
            Idx_s(i,5) = randperm(size(D,1),1);
        end
        
    end

D_val_s = D(Idx_s);
co_sum_s = sum(D_val_s,2);
co_sumsum_s = sum(co_sum_s);
P_s = co_sum_s ./ co_sumsum_s;
H_s = -(P_s.*log2(P_s)); 
H_sum_s(j,:)=sum(H_s); % shuffle entropy
clear D_val_s co_sum_s co_sumsum_s P_s H_s Idx_s
end
entropy_s_sort=sort(H_sum_s);
entropy_s_inf95=entropy_s_sort(1e4*0.05+1);

% ＜TEST＞compare real entropy and surrogate
% The smaller the entropy value, the less dispersed
if H_sum < entropy_s_inf95
    disp('SIGNIFICANT!!!');
else
    disp('Not significant..');
end

surrogate_p=entropy_s_sort;

Pv=(([1:1e4]/1e4)');
surrogate_p(:,2)=Pv;

[IDX]=knnsearch(surrogate_p(:,1),H_sum);
ao=round(surrogate_p(IDX),4);
fao=find(round(surrogate_p(:,1),4)==ao);

mfao=max(fao);
pvalue=Pv(mfao);


%% MDS(fig4A)
% 210623
% 210716 update (after cross-validation)
% load x_test from real_cv
% import data(x_test)
csvfiles=dir('*.csv');%be careful for the PATH!!
numfiles=length(csvfiles);
x_test=cell(1,numfiles);

for k=1:numfiles
    x_test{k}=importdata(csvfiles(k).name);
end
clear k numfiles

x_test = x_test';
x_test = cell2mat(x_test);

% import data(decoded_vm)
% csvfiles=dir('*.csv');%be careful for the PATH!!
% numfiles=length(csvfiles);
% decoded_vm=cell(1,numfiles);
% 
% for k=1:numfiles
%     decoded_vm{k}=importdata(csvfiles(k).name);
% end
% clear k numfiles
% 
decoded_vm = decoded_vm';
decoded_vm = cell2mat(decoded_vm);

% if there is duplicate in x_test, MDS cant be executed. → delete one in that case
for j= 1:size(x_test,1)
    for i=1:size(x_test,1)
        F(i,j) = isequal(x_test(j,:) , x_test(i,:));
    end
end
clear i j

FF = triu(F,1);
FFsum = sum(FF,2);
index_exclude = find(FFsum == 1); %duplicated idx

x_test(index_exclude,:)=[];
decoded_vm(index_exclude,:)=[];

% similar work to shuffled data（(SWnum*10)×2001 matrix）
for j=1:size(index_exclude,1);
    for i=1:9;
        idx_ex(j,i)=index_exclude(j)+(180*i);
    end
end
clear i j
idx_ex = [index_exclude idx_ex];
idx_ex = reshape(idx_ex,[],1);

% for 3cells
RMSE_corr_cv(index_exclude,:)=[];
RMSE_corr_cv_s(index_exclude,:)=[];
RMSE_corr_cv_3c1(index_exclude,:)=[];
RMSE_corr_cv_3c1_s(index_exclude,:)=[];
RMSE_corr_cv_3c2(index_exclude,:)=[];
RMSE_corr_cv_3c2_s(index_exclude,:)=[];
clear RMSE_corr_cv_s_m RMSE_corr_cv_3c1_s_m RMSE_corr_cv_3c2_s_m

RMSE_corr_cv_s_m = mean(RMSE_corr_cv_s,2);
RMSE_corr_cv_3c1_s_m = mean(RMSE_corr_cv_3c1_s,2);
RMSE_corr_cv_3c2_s_m = mean(RMSE_corr_cv_3c2_s,2);

% for 4cells
% RMSE_corr_cv(index_exclude,:)=[];
% RMSE_corr_cv_s(index_exclude,:)=[];
% RMSE_corr_cv_4c1(index_exclude,:)=[];
% RMSE_corr_cv_4c1_s(index_exclude,:)=[];
% RMSE_corr_cv_4c2(index_exclude,:)=[];
% RMSE_corr_cv_4c2_s(index_exclude,:)=[];
% RMSE_corr_cv_4c3(index_exclude,:)=[];
% RMSE_corr_cv_4c3_s(index_exclude,:)=[];
% clear RMSE_corr_cv_s_m RMSE_corr_cv_4c1_s_m RMSE_corr_cv_4c2_s_m RMSE_corr_cv_4c3_s_m
% 
% RMSE_corr_cv_s_m = mean(RMSE_corr_cv_s,2);
% RMSE_corr_cv_4c1_s_m = mean(RMSE_corr_cv_4c1_s,2);
% RMSE_corr_cv_4c2_s_m = mean(RMSE_corr_cv_4c2_s,2);
% RMSE_corr_cv_4c3_s_m = mean(RMSE_corr_cv_4c3_s,2);

% for 5cells
% RMSE_corr_cv(index_exclude,:)=[];
% RMSE_corr_cv_s(index_exclude,:)=[];
% RMSE_corr_cv_5c1(index_exclude,:)=[];
% RMSE_corr_cv_5c1_s(index_exclude,:)=[];
% RMSE_corr_cv_5c2(index_exclude,:)=[];
% RMSE_corr_cv_5c2_s(index_exclude,:)=[];
% RMSE_corr_cv_5c3(index_exclude,:)=[];
% RMSE_corr_cv_5c3_s(index_exclude,:)=[];
% RMSE_corr_cv_5c4(index_exclude,:)=[];
% RMSE_corr_cv_5c4_s(index_exclude,:)=[];
% clear RMSE_corr_cv_s_m RMSE_corr_cv_5c1_s_m RMSE_corr_cv_5c2_s_m RMSE_corr_cv_5c3_s_m RMSE_corr_cv_5c4_s_m
% 
% RMSE_corr_cv_s_m = mean(RMSE_corr_cv_s,2);
% RMSE_corr_cv_5c1_s_m = mean(RMSE_corr_cv_5c1_s,2);
% RMSE_corr_cv_5c2_s_m = mean(RMSE_corr_cv_5c2_s,2);
% RMSE_corr_cv_5c3_s_m = mean(RMSE_corr_cv_5c3_s,2);
% RMSE_corr_cv_5c4_s_m = mean(RMSE_corr_cv_5c4_s,2);


% creating input arguments required for MDS execution（every RMSE in each SW）
for i=1:size(x_test,1)
    for j=1:size(x_test,1)
        RMSE_xtest(i,j) = sqrt(immse(x_test(i,:),x_test(j,:)));
    end
end
clear i j

Y = mdscale(RMSE_xtest,2,'criterion','metricstress');
%%

% 3cells
shuffle_mean = [mean(RMSE_corr_cv_3c1_s(:,1:100),2) mean(RMSE_corr_cv_3c1_s(:,101:200),2) mean(RMSE_corr_cv_3c1_s(:,201:300),2)];

% % 4cells
% shuffle_mean = [mean(RMSE_corr_cv_4c1_s(:,1:100),2) mean(RMSE_corr_cv_4c1_s(:,101:200),2) mean(RMSE_corr_cv_4c1_s(:,201:300),2) mean(RMSE_corr_cv_4c1_s(:,301:400),2)];
% 
% % 5cells
% shuffle_mean = [mean(RMSE_corr_cv_5c1_s(:,1:100),2) mean(RMSE_corr_cv_5c1_s(:,101:200),2) mean(RMSE_corr_cv_5c1_s(:,201:300),2) mean(RMSE_corr_cv_5c1_s(:,301:400),2) mean(RMSE_corr_cv_5c1_s(:,401:500),2)];

%★(S1-R1)/S1
for i = 1:size(shuffle_mean,1)
    for j= 1:size(shuffle_mean,2)
        R1_real(i,j) = shuffle_mean(i,j) - RMSE_corr_cv_3c1(i,j); 
        R1_real_rate(i,j) = R1_real(i,j)/shuffle_mean(i,j);
    end
end

R1_real_rate_1=R1_real_rate(:,1);
R1_real_rate_2=R1_real_rate(:,2);
R1_real_rate_3=R1_real_rate(:,3);
% R1_real_rate_4=R1_real_rate(:,4);
% R1_real_rate_5=R1_real_rate(:,5);


[~,ind1] = sort(R1_real_rate_1);
[~,ind2] = sort(R1_real_rate_2);
[~,ind3] = sort(R1_real_rate_3);
% [~,ind4] = sort(R1_real_rate_4);
% [~,ind5] = sort(R1_real_rate_5);


% entropy -- (S1-R1)/S1 of each slice --
% 210624

% Take log2 to (S1-R1)/S1 for calculation of entropy.
% The values of (S1-R1)/S1 indicated near 0 most in cases so the probability tends to be quite small. 
R1_real_rate_1_log = real(log2(R1_real_rate(:,1)));
R1_real_rate_2_log = real(log2(R1_real_rate(:,2)));
R1_real_rate_3_log = real(log2(R1_real_rate(:,3)));
% R1_real_rate_4_log = real(log2(R1_real_rate(:,4)));
% R1_real_rate_5_log = real(log2(R1_real_rate(:,5)));

[Idx,~] = knnsearch(Y,Y,'K',5);
clear ch_rate1 ch_rate1_sum ch_rate1_sum ch_rate1_sumsum P1 H1 H1_sum
ch_rate1 = R1_real_rate_1_log(Idx);
ch_rate1_sum = sum(ch_rate1,2);
ch_rate1_sumsum = sum(ch_rate1_sum);
P1 = ch_rate1_sum ./ ch_rate1_sumsum;
H1 = -(P1.*log2(P1));
H1 = real(H1);
H1_sum = sum(H1); % real entropy

clear ch_rate2 ch_rate2_sum ch_rate2_sum ch_rate2_sumsum P2 H2 H2_sum
ch_rate2 = R1_real_rate_2_log(Idx);
ch_rate2_sum = sum(ch_rate2,2);
ch_rate2_sumsum = sum(ch_rate2_sum);
P2 = ch_rate2_sum ./ ch_rate2_sumsum;
H2 = -(P2.*log2(P2));
H2 = real(H2);
H2_sum = sum(H2); % real entropy

clear ch_rate3 ch_rate3_sum ch_rate3_sum ch_rate3_sumsum P3 H3 H3_sum
ch_rate3 = R1_real_rate_3_log(Idx);
ch_rate3_sum = sum(ch_rate3,2);
ch_rate3_sumsum = sum(ch_rate3_sum);
P3 = ch_rate3_sum ./ ch_rate3_sumsum;
H3 = -(P3.*log2(P3));
H3 = real(H3);
H3_sum = sum(H3); % real entropy
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

% shuffle Idx（N×5）
clear I
I = [1:size(x_test,1)]';
IDX={};
for j=1:1e6
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
    idx = num2cell(Idx_s,1);
    IDX=[IDX;idx];

ch_rate1_s = R1_real_rate_1_log(Idx_s);
ch_rate1_sum_s = sum(ch_rate1_s,2);
ch_rate1_sumsum_s = sum(ch_rate1_sum_s);
P1_s = ch_rate1_sum_s ./ ch_rate1_sumsum_s;
H1_s = -(P1_s.*log2(P1_s));
H1_s = real(H1_s);
H1_sum_s(j,:)=sum(H1_s); % shuffle entropy

ch_rate2_s = R1_real_rate_2_log(Idx_s);
ch_rate2_sum_s = sum(ch_rate2_s,2);
ch_rate2_sumsum_s = sum(ch_rate2_sum_s);
P2_s = ch_rate2_sum_s ./ ch_rate2_sumsum_s;
H2_s = -(P2_s.*log2(P2_s));
H2_s = real(H2_s);
H2_sum_s(j,:)=sum(H2_s); % shuffle entropy

ch_rate3_s = R1_real_rate_3_log(Idx_s);
ch_rate3_sum_s = sum(ch_rate3_s,2);
ch_rate3_sumsum_s = sum(ch_rate3_sum_s);
P3_s = ch_rate3_sum_s ./ ch_rate3_sumsum_s;
H3_s = -(P3_s.*log2(P3_s));
H3_s = real(H3_s);
H3_sum_s(j,:)=sum(H3_s); % shuffle entropy
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
entropy_s_sort2=sort(H2_sum_s);
entropy_s_sort3=sort(H3_sum_s);
% entropy_s_sort4=sort(H4_sum_s);
% entropy_s_sort5=sort(H5_sum_s);

clear entropy_s_inf95_1 entropy_s_inf95_2 entropy_s_inf95_3 entropy_s_inf95_4 entropy_s_sup95_1 entropy_s_sup95_2 entropy_s_sup95_3
entropy_s_inf95_1=entropy_s_sort1(1e4*0.05+1);
entropy_s_inf95_2=entropy_s_sort2(1e4*0.05+1);
entropy_s_inf95_3=entropy_s_sort3(1e4*0.05+1);
% entropy_s_inf95_4=entropy_s_sort4(1e4*0.05+1);
% entropy_s_inf95_5=entropy_s_sort5(1e4*0.05+1);


% ＜TEST＞compare real entropy and surrogate
% The smaller the entropy value, the less dispersed

if H1_sum < entropy_s_inf95_1
    disp('SIGNIFICANT!!!');
else
    disp('Not significant...');
end
if H2_sum < entropy_s_inf95_2
    disp('SIGNIFICANT!!!');
else
    disp('Not significant...');
end
if H3_sum < entropy_s_inf95_3
    disp('SIGNIFICANT!!!');
else
    disp('Not significant...');
end
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
surrogate_p2=entropy_s_sort2;
surrogate_p3=entropy_s_sort3;
% surrogate_p4=entropy_s_sort4;
% surrogate_p5=entropy_s_sort5;

clear Pv
Pv=(([1:1e6]/1e6)');
surrogate_p1(:,2)=Pv;
surrogate_p2(:,2)=Pv;
surrogate_p3(:,2)=Pv;
% surrogate_p4(:,2)=Pv;
% surrogate_p5(:,2)=Pv;

clear IDX1 IDX2 IDX3 IDX4 IDX5
[IDX1]=knnsearch(surrogate_p1(:,1),H1_sum);
[IDX2]=knnsearch(surrogate_p2(:,1),H2_sum);
[IDX3]=knnsearch(surrogate_p3(:,1),H3_sum);
% [IDX4]=knnsearch(surrogate_p4(:,1),H4_sum);
% [IDX5]=knnsearch(surrogate_p5(:,1),H5_sum);

clear ao1 ao2 ao3 ao4 ao5 fao1 fao2 fao3 fao4 fao5 mfao1 mfao2 mfao3 mfao4 mfao5 pvalue_all
ao1=round(surrogate_p1(IDX1),4);
ao2=round(surrogate_p2(IDX2),4);
ao3=round(surrogate_p3(IDX3),4);
% ao4=round(surrogate_p4(IDX4),4);
% ao5=round(surrogate_p5(IDX5),4);

fao1=find(round(surrogate_p1(:,1),4)==ao1);
fao2=find(round(surrogate_p2(:,1),4)==ao2);
fao3=find(round(surrogate_p3(:,1),4)==ao3);
% fao4=find(round(surrogate_p4(:,1),4)==ao4);
% fao5=find(round(surrogate_p5(:,1),4)==ao5);


mfao1=max(fao1);
mfao2=max(fao2);
mfao3=max(fao3);
% mfao4=max(fao4);
% mfao5=max(fao5);


pvalue_all=[Pv(mfao1) Pv(mfao2) Pv(mfao3)]'% Pv(mfao4) Pv(mfao5)



