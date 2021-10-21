%% ripple_finalization2.m
%% 190108 update, Nobuyoshi Matsumoto, Univ Tokyo, 2019
% ripple_finalization.mだとonsetとoffsetが微妙なので2SDで取ることにした。
% ripple_finalization.mの後に使う。
% 保存されるファイルのうち、ファイル名がfinal2となっているものが、このプログラムで求めた変数や図

%% USAGE
% EXAMPLE: ripple_finalization2('D:\data\nm0006\171208');
%
% Input:
% - basepath: ripple or noiseを目で振り分けたExcelを保存しているフォルダ
%   basepathに図や計算結果が保存される。
%   basepathにmatファイルが一つしかない場合はそれを読み込む。
%   basepathにatfファイルが複数ある場合は、ユーザーインターフェイスから選択して読み込む。
% Output:
% datasaveoptionを1にしておくと、計算結果（AllResults）と設定（AllSettings）がbasepathに保存される。
% figsaveoptionを1にしておくと、図（.figと.tif）がbasepathに保存される。
% deleteoldresultsを1にしておくと、ripple_candidates_detection.mで計算された古い結果を消去する。

function ripple_finalization2(basepath)
tic; cd(basepath); clearvars -except basepath

%% Set patameters
datasaveoption = 1; % 1の場合、変数（計算結果）を保存する。
figsaveoption = 0; % 1の場合、図（.fig）を保存する。今回の場合、.tifより重い！
tifsaveoption = 1; % 1の場合、tiffの図を保存する。.figファイルよりは軽い。
deleteoldresults = 0; % 1の場合、ripple_candidates_detection.mで計算された古い結果を消去する。

ThrONOFF = 2; % オンセットとオフセットを決めるための閾値
FminSGamma = 20; % スローガンマ帯域の下限周波数
FmaxSGamma = 40; % スローガンマ帯域の上限周波数

%% Load MAT file
matfile = dir(fullfile(basepath,'*_final.mat'));
if length(matfile) < 1; error('No mat file detected in the folder!'); return;
elseif length(matfile) == 1; matname = matfile.name;
else
    [matname,fpath]=uigetfile('*.mat','Select MAT file of ripple CANDIDATES');
    if isequal([matname fpath],[0 0])
        display('User canceled.');
        return;
    end
end
disp(strcat('Loading:',matname));
load(matname);

ARfs       = AllResults.fs;
ARdataname = AllResults.dataname;

NewDataName = strcat(ARdataname(1:end-10),'_final2.mat');
NewTifName  = strcat(ARdataname(1:end-10),'_final2.tif');
NewFigName  = strcat(ARdataname(1:end-10),'_final2.fig');

%% Load XLSX file
% -- Not needed

%% Derive wave information (in the almost same way as a part of ripple_candidates_detection.m)
% -- Import ATF file
atffile = dir(fullfile(basepath,'*.atf'));
if length(atffile) < 1; error('No atf file detected in the folder!'); return;
elseif length(atffile) == 1; atfname = atffile.name;
else
    [atfname,fpath]=uigetfile('*.atf','Select ATF file');
    if isequal([atfname fpath],[0 0])
        display('User canceled.');
        return;
    end
end
disp(strcat('Loading...',atfname));
[~,~,~,d]=import_atf(atfname);
disp(strcat('Loading done!',atfname));
ARtime=d(:,1);
ARlfp0=d(:,AllSettings.ii); %LFPが含まれる列を指定（LFPの波形）
fs=ARfs; % sampling frequency
clear d

% -- Filter LFP traces
[B,A]=butter(3,[AllSettings.Fmin/(fs*0.5) AllSettings.Fmax/(fs*0.5)]);
ARlfpf=filtfilt(B,A,ARlfp0);

% -- Calculate RMS
rmsdurainidx=round(AllSettings.rms_duration*fs/1000);
forenlarge_y=zeros(rmsdurainidx,1);
forenlarge_y(:)=mean(ARlfpf); % 生波形にフィルターをかけた波形の電位の平均値を縦にrmsdurationidx個並べた。
temp_y=[forenlarge_y;ARlfpf;forenlarge_y]; % 両端を補正するために、最初と最後にforenlarge_yをくっつけた。
temp_y=temp_y-mean(temp_y); % baselineを0にした。temp_yからmeanを引いても、先にtemp_yを平均して、前後に0の列をくっつけても結果は同じ。
ARrms_y=zeros(size(ARlfpf));

for i=1:size(ARlfpf(:,1))   %二乗平均平方根，RMS
   ARrms_y(i)=sqrt(sum(temp_y(i:i+rmsdurainidx-1).^2)/rmsdurainidx);   
end
clear i

Zripplerms=(ARrms_y-mean(ARrms_y)*ones(size(ARrms_y)))./(std(ARrms_y)*ones(size(ARrms_y)));

%% Detect ripple peaks: find max RMS
% peakを求め直す必要はない。

%% Detect ripple onsets and offsets
ONOFFidx=zeros(size(ARlfp0)); % リップルの候補はここに1を入れていく（最初はzerosで全て0にしておく）
% ARlfp0が列ベクトルなら、rippleidxも列ベクトルになる。
ONOFFidx(Zripplerms>=ThrONOFF)=1; % thresholdを超えた部分に1を入れる。
DiffONOFFidx=[ONOFFidx(1);diff(ONOFFidx)];
ripplefinalindex2=AllResults.ripplefinalindex; % onsetとoffsetを求め直し、かつ5列に、ピークが何SDかを代入した
% -- Re-find onsets and offsets --
ON  = find(DiffONOFFidx==1);
OFF = find(DiffONOFFidx==-1)-1;
allripplenum=size(ripplefinalindex2,1);
for qq=1:allripplenum
    ripplefinalindex2(qq,2) =...
        ON (max(find( ON<AllResults.ripplefinalindex(qq,2)))); % Find onsets
    ripplefinalindex2(qq,4) =...
        OFF(min(find(OFF>AllResults.ripplefinalindex(qq,4)))); % Find offsets
end
clear qq

ripplefinaltiming2=[ripplefinalindex2(:,1),ARtime(ripplefinalindex2(:,2:4))];

%% Calculate other parameters
RippleInfo2=[ripplefinalindex2,ripplefinaltiming2,...
             ripplefinaltiming2(:,4)-ripplefinaltiming2(:,2),...
             ((ripplefinalindex2(:,3)-ripplefinalindex2(:,2))-(ripplefinalindex2(:,4)-ripplefinalindex2(:,3)))./...
             (ripplefinalindex2(:,4)-ripplefinalindex2(:,2)),...
             zeros(allripplenum,3),...
             Zripplerms(ripplefinalindex2(:,3))]; % この時点では14列にしておく。
%  1列　リップル番号
%  2列　オンセットのインデックス
%  3列　ピークのインデックス
%  4列　オフセットのインデックス
%  5列　リップル番号
%  6列　オンセットの時刻（秒）
%  7列　ピークの時刻（秒）
%  8列　オフセットの時刻（秒）
%  9列　持続時間（秒）
% 10列　Symmetry index (=(a-b)/(a+b): a ピークとオンセットの距離、b オフセットとピークの距離）
% 11列　リップル帯域のウェーブレットパワー（ウェーブレットスペクトルの強度の和を時間で規格化、Fmin～Fmax Hz）
%   　　　…取ってくる時間帯が短すぎてFFTだとスペクトルが粗すぎるので、ウェーブレットにした
% 12列　スローガンマ帯域のFFTパワー（ウェーブレットスペクトルの強度の和を時間で規格化、20～40 Hz）
% 13列　リップル帯域のピーク周波数（ウェーブレットスペクトルのリップル帯域において、PSDが最大になる周波数）
% 14列　リップルのピークが何SDか

f_axis=[1:.01:2.5];
r1 =find(f_axis>log10(AllSettings.Fmin),1,'first'); % リップルの最低周波数の常用対数のf_axisにおける「インデックス」
r2 =find(f_axis<log10(AllSettings.Fmax),1, 'last');
sg1=find(f_axis>log10(FminSGamma),1,'first');
sg2=find(f_axis>log10(FmaxSGamma),1, 'last');

for qq=1:allripplenum
    G=ARlfp0(ripplefinalindex2(qq,2):ripplefinalindex2(qq,4));
    [coef, abscoef, time]=wavelet_nm_log_nofig(G,0,length(G)/fs,fs,f_axis,'cmor1.5-2');
    RippleInfo2(qq,11)=sum(sum(abscoef( r1: r2,:)))/(length(G)/fs);
    RippleInfo2(qq,12)=sum(sum(abscoef(sg1:sg2,:)))/(length(G)/fs);
    Sp=sum(abscoef,2); SpRip=Sp(r1:r2); IdxMaxRipplePower=find(SpRip==max(SpRip));
    RippleInfo2(qq,13)=10^f_axis((r1-1)+IdxMaxRipplePower);
end


%% Draw figures (and save if required)
% 描画 
SS = get(0, 'ScreenSize'); % スクリーンサイズを取得
figure('Position',SS);
subplot(211); % raw (non-filtered) trace
plot(ARtime,ARlfp0,'k'); hold on; axis tight
scatter(AllResults.ripplefinaltiming(:,2),ARlfp0(AllResults.ripplefinalindex(:,2)),60,'>','MarkerEdgeColor','b');
scatter(AllResults.ripplefinaltiming(:,4),ARlfp0(AllResults.ripplefinalindex(:,4)),60,'<','MarkerEdgeColor','b');
scatter(ripplefinaltiming2(:,2),ARlfp0(ripplefinalindex2(:,2)),60,'>','MarkerEdgeColor','g');
scatter(ripplefinaltiming2(:,4),ARlfp0(ripplefinalindex2(:,4)),60,'<','MarkerEdgeColor','g');
scatter(AllResults.ripplefinaltiming(:,3),ARlfp0(AllResults.ripplefinalindex(:,3)),60,'o','MarkerEdgeColor','r');
for i=1:allripplenum
    plot(ARtime(ripplefinalindex2(i,2):ripplefinalindex2(i,4)),...
         ARlfp0(ripplefinalindex2(i,2):ripplefinalindex2(i,4)),'c');
end
for i=1:allripplenum
    plot(ARtime(AllResults.ripplefinalindex(i,2):AllResults.ripplefinalindex(i,4)),...
         ARlfp0(AllResults.ripplefinalindex(i,2):AllResults.ripplefinalindex(i,4)),'m');
end
title(['Raw LFP trace after visual inspection (Data: ',atfname,'; Data sampling: ',...
    num2str(fs),' Hz; Ripple detection threshold: ',num2str(AllSettings.ThrSD),'SD of RMS of filtered trace (',...
    num2str(AllSettings.Fmin),'-',num2str(AllSettings.Fmax),' Hz))'],...
    'interpreter','none');
% 'interpreter','none'
% アンダーバーをそのままタイトルに表示できるようになる
%（これを入れないと、アンダーバーの次の文字が下付きになる）
% 最後のReferenceも参考に。
xlabel('Time (sec)');
hold off;

subplot(212); % filtered trace
plot(ARtime,ARlfpf,'k'); hold on; axis tight
scatter(AllResults.ripplefinaltiming(:,2),ARlfpf(AllResults.ripplefinalindex(:,2)),60,'>','MarkerEdgeColor','b');
scatter(AllResults.ripplefinaltiming(:,4),ARlfpf(AllResults.ripplefinalindex(:,4)),60,'<','MarkerEdgeColor','b');
scatter(ripplefinaltiming2(:,2),ARlfpf(ripplefinalindex2(:,2)),60,'>','MarkerEdgeColor','g');
scatter(ripplefinaltiming2(:,4),ARlfpf(ripplefinalindex2(:,4)),60,'<','MarkerEdgeColor','g');
scatter(AllResults.ripplefinaltiming(:,3),ARlfpf(AllResults.ripplefinalindex(:,3)),60,'o','MarkerEdgeColor','r');
for i=1:allripplenum
    plot(ARtime(ripplefinalindex2(i,2):ripplefinalindex2(i,4)),...
         ARlfpf(ripplefinalindex2(i,2):ripplefinalindex2(i,4)),'c');
end
for i=1:allripplenum
    plot(ARtime(AllResults.ripplefinalindex(i,2):AllResults.ripplefinalindex(i,4)),...
         ARlfpf(AllResults.ripplefinalindex(i,2):AllResults.ripplefinalindex(i,4)),'m');
end
title(['Filtered LFP trace after visual inspection (Data: ',atfname,'; Data sampling: ',...
    num2str(fs),' Hz; Threshold: ',num2str(AllSettings.ThrSD),'SD, Bandpass: ',...
    num2str(AllSettings.Fmin),'-',num2str(AllSettings.Fmax),' Hz)'],...
    'interpreter','none');
xlabel('Time (sec)');
hold off;

if figsaveoption
    saveas(gcf,NewFigName,'fig');
end

if tifsaveoption
    saveas(gcf,NewTifName,'tiffn');
end

%% Delete old data and figures if necessary
if deleteoldresults % 最初の設定で、deleteoldresultsを1にしておくと、古い変数や図を消去できる。
    delete([ARdataname(1:end-4),'.mat']);
    delete([ARdataname(1:end-4),'.tif']);
    delete([ARdataname(1:end-4),'.fig']);
end

%% Save results and settings
% 変数や設定を保存する
cd(basepath);
if datasaveoption
    AllResults.ripplefinaltiming2 = ripplefinaltiming2;
    AllResults.ripplefinalindex2  = ripplefinalindex2;
    AllResults.fs         = ARfs;
    AllResults.dataname   = NewDataName;
    AllResults.allripplenum   = allripplenum;
    AllResults.RippleInfo2 = RippleInfo2;
    
    AllSettings.ThrONOFF = ThrONOFF;
    AllSettings.FminSGamma = FminSGamma;
    AllSettings.FmaxSGamma = FmaxSGamma;

    clearvars -except AllResults AllSettings
    save(AllResults.dataname);
end

toc;

%% Reference
% タイトル等のTextオブジェクトにアンダーバー「_」やハット「^」を表示したいのですが、
% それぞれの記号の直後の文字が下付?き文字と上付き文字になってしまいます。
% https://jp.mathworks.com/matlabcentral/answers/98127-text-_
%
% findをfirstやlastと併用する
% https://jp.mathworks.com/help/matlab/ref/find.html
% k = find(X,n,direction) は、direction が 'last' の場合、
% X 内の非ゼロ要素に対応するインデックスを最後から数えて n 個検出します。
% direction の既定値は 'first' であり、非ゼロ要素に対応するインデックスを最初から数えて n 個検出します。
