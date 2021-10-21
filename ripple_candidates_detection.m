%% ripple_candidates_detection.m
%% 18803 update, Nobuyoshi Matsumoto, OCU, 2018

%% USAGE
% EXAMPLE: ripple_detection_0('D:\data\nm0006\171208');
%
% Input:
% - basepath: atf保存しているフォルダ。例：'D:\data\nm0006\171208'
%   basepathに図や計算結果（およびExcel）が保存される。
%   basepathにatfファイルが一つしかない場合はそれを読み込む。
%   basepathにatfファイルが複数ある場合は、ユーザーインターフェイスから選択して読み込む。
% Output:
% datasaveoptionを1にしておくと、計算結果（AllCandidates）と設定（AllSettings）がbasepathに保存される。
% figsaveoptionを1にしておくと、図（.figと.tif）がbasepathに保存される。
% xlssaveoptionを1にしておくと、リップル／ノイズの振り分け用のxlsがbasepathに保存される。

function ripple_candidates_detection(basepath)
tic;
cd(basepath);
clearvars -except basepath

%% Set patameters
ii=3; % importしたatfのデータの何列にLFPのデータが入っているか。defaultは3??

Fmax = 249; % フィルタ周波数最大値
Fmin = 120; % フィルタ周波数最小値
ThrSD = 3; % RMSの閾値（ベースラインから何SDか）)

rms_duration = 3; % in msec　　rmsを計算するためのbinのduration
min_duration = 10; % in msec   これを下回る期間をイベントと判定しない。
min_isw_period = 5; % in msec ripple間がこれより短ければmergeする。

datasaveoption = 1; % 1の場合、変数（計算結果）を保存する。
figsaveoption = 0; % 1の場合、図（拡張子 .fig・・・重い！）を保存する。この段階では保存しなくても良い（リップルが確定ではないし重いので）。
tifsaveoption = 1; % 1の場合、tiffの図を保存する。.figファイルよりは軽い。
xlsmakeoption = 1; % 1の場合、Excelファイルを作成する。

%% Select and import ATF file
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
[a,b,c,d]=import_atf(atfname); clear a b c
disp(strcat('Loading:',atfname));

time=d(:,1);
lfp0=d(:,ii); %LFPが含まれる列を指定（LFPの波形）
fs=round(1/(d(2,1)-d(1,1))); % sampling frequency
clear d

%% Set file names
filenamebase  = strcat(atfname(1:end-4),'_',num2str(Fmin),'-',num2str(Fmax),'Hz_',num2str(ThrSD),'SD');
xlsname  = [filenamebase,'.xlsx'];
dataname = [filenamebase,'.mat'];
tifname  = [filenamebase,'.tif'];
figname  = [filenamebase,'.fig'];

%% Filter LFP traces
[B,A]=butter(3,[Fmin/(fs*0.5) Fmax/(fs*0.5)]);
y=filtfilt(B,A,lfp0); % filter後のLFPの波形（あとで、yはlfpfと名前を変える）

%% Detect ripple candidates (onset & offset): calculate RMS, threshold, zscore
rmsdurainidx=round(rms_duration*fs/1000); % rms_durationはインデックスに換算するといくつ分になるか。（3や4などの数字ひとつのスカラーが出てくる）
forenlarge_y=zeros(rmsdurainidx,1);
forenlarge_y(:)=mean(y);
temp_y=[forenlarge_y;y;forenlarge_y];
temp_y=temp_y-mean(temp_y); % baselineを0にした。
rms_y=zeros(size(y));

for i=1:size(y(:,1))   %二乗平均平方根，RMS
   rms_y(i)=sqrt(sum(temp_y(i:i+rmsdurainidx-1).^2)/rmsdurainidx);   
end

mean_rms_y=zeros(size(rms_y)); mean_rms_y(:)=mean(rms_y);
std_rms_y=zeros(size(rms_y)); std_rms_y(:)=std(rms_y);

zrippledetection=(rms_y-mean_rms_y)./std_rms_y;

rippleidx=zeros(size(lfp0)); % リップルの候補はここに1を入れていく（最初はzerosで全て0にしておく）
% lfp0が列ベクトルなら、rippleidxも列ベクトルになる。

rippleidx(zrippledetection>=ThrSD)=1; % thresholdを超えた部分に1を入れる。

%% Remove ripple candidates if each duration is very short
% リップルのdurationが短すぎるものは0にする（min_durationを使う）
diffrippleidx=[0;diff(rippleidx)]; % 差をとっただけ。列ベクトルになるはず。
onsetidx =  find(diffrippleidx== 1); % onsetの候補
offsetidx= (find(diffrippleidx==-1))-1; % offsetの候補

numon = numel(onsetidx);
numoff= numel(offsetidx);
% ↑numonとnumoffは同じ数であって、かつ、onsetidxは常にoffsetidxよりも小さくなければならないので、
% 以下のようにする。
if numon > numoff
    onsetidx(end)=[];
elseif numon < numoff
    offsetidx(1)=[];
end

rippleinfo = [[1:numel(onsetidx)]', onsetidx, offsetidx,...
              onsetidx/fs, offsetidx/fs];
% 1列にripple番号、2列にripple onsetのindex、3列にripple offsetのindex、
% 4列にripple onsetの時刻（秒）、5列にripple offsetの時刻（秒）。

for i=1:size(rippleinfo,1) % 短すぎるリップルは削除（まず0で置き換えて、for文が終わったら削除する）
    if (rippleinfo(i,5)-rippleinfo(i,4)) <= min_duration/1000 % msec⇒secにすることに注意
        rippleinfo(i,:)=[0 0 0 0 0];
    end
end
i=rippleinfo(:,1); excludeidx=find(i==0);
rippleinfo(excludeidx,:)=[];
rippleinfo(:,1)=[1:size(rippleinfo,1)]'; % ripple番号を振り直した（1,2,3,4...になるようにした）
rippleinfo=[rippleinfo,zeros(size(rippleinfo,1),3)]; % 6列にrippleのduration,7列にpeakのindex, 8列にpeakの時刻を入れる
rippleinfo(:,6)=(rippleinfo(:,5)-rippleinfo(:,4))*1000; % duration、単位は「ミリ秒」

%% Merge ripple candidates if inter-ripple intarvals are very short
% 近すぎるリップルはmergeする
i=1;
while i < size(rippleinfo,1)
    if ~isempty(rippleinfo(i,1))
        if rippleinfo(i+1,2)-rippleinfo(i,3) < (min_isw_period/1000)*fs
            rippleinfo(i,3)=rippleinfo(i+1,3);
            rippleinfo(i,5)=rippleinfo(i+1,5);
            rippleinfo(i+1,:)=[];
        end
    end
i=i+1;
end
rippleinfo(:,1)=[1:size(rippleinfo,1)]';

%% Detect ripple candidates (peak): find max RMS
for i=1:size(rippleinfo,1) % ピークの検出
    eachrms=rms_y(rippleinfo(i,2):rippleinfo(i,3));
    pkidxcand=find(eachrms==max(eachrms),1,'first');
    rippleinfo(i,7)=rippleinfo(i,2)+pkidxcand-1; %-1をしないとピークのインデックスからひとつずれてしまう
    rippleinfo(i,8)=time(rippleinfo(i,7));
end

rippleinfo(:,4)=time(rippleinfo(:,2));
rippleinfo(:,5)=time(rippleinfo(:,3));

%% Draw figures (and save if required)
% 描画 
SS = get(0, 'ScreenSize'); % スクリーンサイズを取得
figure('Position',SS);
subplot(211); % raw (non-filtered) trace
plot(time,lfp0,'k'); hold on; axis tight
scatter(time(rippleinfo(:,2)),lfp0(rippleinfo(:,2)),'>','MarkerEdgeColor','b');
scatter(time(rippleinfo(:,3)),lfp0(rippleinfo(:,3)),'<','MarkerEdgeColor','b');
scatter(time(rippleinfo(:,7)),lfp0(rippleinfo(:,7)),'o','MarkerEdgeColor','r');
for i=1:size(rippleinfo,1)
    plot(time(rippleinfo(i,2):rippleinfo(i,3)),lfp0(rippleinfo(i,2):rippleinfo(i,3)),'m');
end
title(['Raw LFP trace before visual inspection (',num2str(ThrSD),'SD, ',...
    num2str(Fmin),'-',num2str(Fmax),' Hz)']);
xlabel('Time (sec)');
hold off;
subplot(212); % filtered trace
plot(time,y,'k'); hold on; axis tight
scatter(time(rippleinfo(:,2)),y(rippleinfo(:,2)),'>','MarkerEdgeColor','b');
scatter(time(rippleinfo(:,3)),y(rippleinfo(:,3)),'<','MarkerEdgeColor','b');
scatter(time(rippleinfo(:,7)),y(rippleinfo(:,7)),'o','MarkerEdgeColor','r');
for i=1:size(rippleinfo,1)
    plot(time(rippleinfo(i,2):rippleinfo(i,3)),y(rippleinfo(i,2):rippleinfo(i,3)),'m');
end
title(['Filtered LFP trace before visual inspection (',num2str(ThrSD),'SD, ',...
    num2str(Fmin),'-',num2str(Fmax),' Hz)']);
xlabel('Time (sec)');
hold off;

if figsaveoption
    saveas(gcf,figname,'fig');
end

if tifsaveoption
    saveas(gcf,tifname,'tiffn');   
end

%% Make XLS file
% ここで、rippleの振り分け（ノイズ or シグナル）用のxlsファイルを作る。
if xlsmakeoption
    xlsinfo=[rippleinfo(:,1), zeros(size(rippleinfo,1),1),...
             rippleinfo(:,4), rippleinfo(:,8), rippleinfo(:,5)];
         % 1列にripple番号、2列に0 or 1 or 2を書き込んでいく。
         % 3列にripple onsetの時刻（秒）、4列にripple peakの時刻（秒）、
         % 5列にripple offsetの時刻（秒）。
    xlswrite(xlsname,xlsinfo); % xlsxファイルをここで作る。
    % このxlsファイルの2列に、数値を入れていく。
    % 0: ノイズ
    % 1: 単独でリップル
    % 2: つなげるリップルの最初
    % 0: つなげるリップルの最初でも最後でもないところ（これはあってもなくてもいい）
    % 3: つなげるリップルの最後（2と同じ個数になるはず）
end

%% Save candidates
% 変数や設定を保存する
cd(basepath);
lfpf = y; % 覚えやすくするために、filter後のlfpの変数の名前をlfpfとしただけ。filtered LFPという意味でlfpfとした。
if datasaveoption
    AllCandidates.rippleinfo = rippleinfo;
    AllCandidates.fs         = fs;
%     AllCandidates.lfp0       = lfp0;    
%     AllCandidates.lfpf       = lfpf;
%     AllCandidates.time       = time;
%     AllCandidates.rms_y      = rms_y;
    AllCandidates.dataname   = dataname;
    AllCandidates.atfname    = atfname;
        
    AllSettings.Fmax  = Fmax;
    AllSettings.Fmin  = Fmin;
    AllSettings.ThrSD = ThrSD;
    AllSettings.rms_duration = rms_duration;
    AllSettings.min_duration = min_duration;
    AllSettings.min_isw_period = min_isw_period;
    AllSettings.ii    = ii;
    AllSettings.basepath = basepath;    
    
    clearvars -except AllCandidates AllSettings
    save(AllCandidates.dataname);
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
