%% SW_detection_in_vitro_nm5,m %%
%% 20160817 filter_spw_peak_num.mを改変
%% to find local maxima of x (>level). 上に突き出したpeakのpositionとamplitude, 数をかぞえる
%% ピークの閾値とベースラインの補正はmfile上で行う
%% data reduction後のsampling frequencyは1000 Hzを推奨
%% y1: local field potentials (mV)
%% Update: 20170522 %%

function [spwpeakpos, spwpeaktime, spwpeakamp, spwnum, spwfrq, spwonsetpos, spwdur, spwend, y1sw]=SW_detection_in_vitro_ao(y1,t,T,y1_6);
%[spwpeakpos, spwpeaktime, spwpeakamp, spwnum, spwfrq, spwonsetpos, spwdur, spwend, y1sw, x]=ao_SW_detection_in_vitro(y1,t,T,y1_6);
% spwpeakamp (μV)
% spwpeaktime (s)
% spwfrq (Hz)

%%%%%%%%%%%%%
%パラメタ設定

p=0;%baselineを補正する
fs=1/(t(2)-t(1)); % sampling frequency (Hz)
Fmax=30; %フィルタ周波数最大値
Fmin=2; %フィルタ周波数最小値
%%%%%%%%%%%%%
%filterをかけます
%%%%%%%%%%%%%%%%%%%%%%%%%%

[B,A]=butter(3,[Fmin/(fs*0.5) Fmax/(fs*0.5)]);
x=filtfilt(B,A,y1); y1sw=x;

% figure;
% subplot(211);
% plot(t,y1,'k');
% subplot(212);
% plot(t,x,'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%全体のSD値、ベースラインノイズのSD値を算出します。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdv=std(x) %SD
stdv=2*std(x) %2SD
stdv=3*std(x) %3SD
stdv=4*std(x) %4SD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baselineを決める
xmin=25;%plotする範囲を適当にきめる
xmax=35;
ymin=-0.01;
ymax=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1000); plot(t,x,'k');
axis([xmin, xmax, ymin, ymax]);
baseline=ginput(2);
basexmin=baseline(1,1);
basexmax=baseline(2,1);

basewave=x([round(basexmin*fs):round(basexmax*fs)],:);%整数値でなくてはなりません
SD=std(basewave);
ave=mean(basewave);
ave+SD
ave+2*SD
ave+3*SD
ave+4*SD
ave+5*SD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level=input('level(mV) >>\n');%閾値の設定（ｍV）

[spwpeaktime,spwpeakpos]=findpeaks(x-p,'MINPEAKHEIGHT',level,'MINPEAKDISTANCE',round(0.05*fs));
% x=x-p;
% idx=find(x>=level);
% d1=diff(x);
% d1=[0;d1];
% d2=d1(idx);
% d2(find(d2>=0))=1;
% d2(find(d2<0))=2;
% d3=diff(d2);    
% pidx=find(d3>0);
% spwpeakpos=idx(pidx); %インデックスを求めた
spwpeakpos=spwpeakpos(spwpeakpos>0.1*fs); spwpeakpos=spwpeakpos(spwpeakpos<(round(t(end))-0.1)*fs);
% ↑nobumossyの中で、peakの0.1秒前から0.1秒後まででtraceを重ねるから、最初と最後の方のピークは不要。
spwpeakamp=x(spwpeakpos)*1000; % spwpeakampはμV単位にした
% spwpeaktime=spwpeakpos/fs;
spwnum=size(spwpeakpos,1);
spwfrq=size(spwpeakpos,1)/round(t(end));
spwpeaktime=t(spwpeakpos);


spwonsetpos=zeros(size(spwpeakpos)); %SPW onsetのindex

for j=1:numel(spwpeakpos); % ピークから50ms前までのfiltered traceにおいて、baseline+sdを超えた点をonsetとする
    AA=[round(spwpeakpos(j)-0.1*fs):round(spwpeakpos(j))]; %ピークから50ms前までのインデックス（50msにしたのは適当だが、この範囲内でonsetを見つけられそう）
    B=find(y1sw(AA)>(ave+SD));
    if size(B,2)>1 % 行ベクトルを列ベクトルにする
        B=B';
    end
    dB=[2;diff(B)]; % 最初の数字を1より大きい数にしておけば、if isequal(dB,ones(size(dB))==1とかを入れなくてよい
    c=max(find(dB>1)); % ピークに一番近い側を取りたいので、maxを入れた
    BB=B(c);
    CC=AA(BB); % これがfiltered traceにおいてSDを超えた初めての時刻のindex
    
    % ピークの直前のトラフも含めてave+sdを超えてしまった場合は仕方ないので、ピークの直前のトラフからSDだけ超えた点をとってくる
    yo=-y1sw(CC:AA(end)); % 波形を反転させてトラフを見つける
    [~,postr]=findpeaks(yo,'MINPEAKDISTANCE',5);% round(0.005*fs));→5ms（roundではうまくいかないので、数値を直接入れる） % trough：トラフ % トラフの前後10msで良さそうだがとりあえず5msにしてみる
    if isempty(postr)==0 % ←ピーク直前のトラフを含めave+sdを超えた場合
        cc=max(postr);
        AAasc=AA(cc:end); % ccより先は増加するようになっている
        BB=min(find(y1sw(AAasc)>(y1sw(AAasc(1))+SD)));
        if isempty(BB)==1
            CC=AA(cc);
        else
            CC=AAasc(BB);
        end
    end
    spwonsetpos(j,1)=CC;
end

spwend=zeros(size(spwpeakpos));
for j=1:numel(spwpeakpos); % ピークから50ms後までのfiltered traceにおいて、peakposからbaseline+SDに戻った点をend
    A2=[round(spwpeakpos(j)):round(spwpeakpos(j)+0.1*fs)]; %ピークから100ms後までのインデックス（100msにしたのは適当だが、この範囲内でendを見つけられそう）
    B2=find(y1sw(A2)<=(ave+SD));
    c2=min(B2);% ピークに一番近い側を取りたいので、minを入れた。つまりピークに最も近いave+SDを下回るindex。%ここまでおｋ
    BB2=A2(c2); %aveを下回るindexの値
%     CC2=A2(BB2); % これがfiltered traceにおいてSDを超えた初めての時刻のindex
%   spwend(j,1)=BB2;
end

spwdur=zeros(size(spwpeakpos));
for k=1:numel(spwpeakpos)
spwdur(k,1)=spwend(k,1)-spwonsetpos(k,1);
end
spwdur=spwdur/fs; %時刻（s）
clear k

%以下のfigureを確認用に描図する際は17－26行を実行し、必要な変数(x)をだす。
figure;
% subplot(211);
plot(t,x,'k'); xlim([0 t(end)]); % axis tight;
hold on;
plot(t(spwpeakpos),x(spwpeakpos),'ro', t(spwonsetpos),x(spwonsetpos),'c^'); hold on;%hold off;% % ,t(spwend),x(spwend),'go'
%ylabel('Voltage (mV)'); title('Filtered LFP traces ('num2str(length(spwpeakpos)) 'SWs (threshold = ' num2str(level) 'mV), 2-30 Hz');
% subplot(212);
plot(t,y1-0.25,'k'); xlim([0 t(end)]); hold on;% axis tight; 
plot(t(spwpeakpos),y1(spwpeakpos)-0.25,'ro', t(spwonsetpos),y1(spwonsetpos)-0.25,'c^'); hold on;%,t(spwend),y1(spwend)-0.25,'go'
plot(T,y1_6-0.5,'k'); xlim([0 T(end)]); hold on;% axis tight; 
plot(T(spwpeakpos*20),y1_6(spwpeakpos*20)-0.5,'ro', T(spwonsetpos*20),y1_6(spwonsetpos*20)-0.5,'c^'); hold off; %,T(spwend*20),y1_6(spwend*20)-0.5,'go'
xlabel('Time (s)'); ylabel('Voltage (mV)'); title('Raw LFP traces');

clear Fmax Fmin p A B

% figure;
% % subplot(211);
% plot(t,x,'k'); xlim([0 t(end)]); % axis tight;
% hold on;
% plot(t(spwpeakpos),x(spwpeakpos),'ro', t(spwonsetpos),x(spwonsetpos),'c^',t(spwend),x(spwend),'go'); hold on;%hold off;% % 
% %ylabel('Voltage (mV)'); title('Filtered LFP traces ('num2str(length(spwpeakpos)) 'SWs (threshold = ' num2str(level) 'mV), 2-30 Hz');
% % subplot(212);
% plot(t,y1-0.25,'k'); xlim([0 t(end)]); hold on;% axis tight; 
% plot(t(spwpeakpos),y1(spwpeakpos)-0.25,'ro', t(spwonsetpos),y1(spwonsetpos)-0.25,'c^',t(spwend),y1(spwend)-0.25,'go'); hold on;%
% plot(T,y1_6-0.5,'k'); xlim([0 T(end)]); hold on;% axis tight; 
% plot(T(spwpeakpos*20),y1_6(spwpeakpos*20)-0.5,'ro', T(spwonsetpos*20),y1_6(spwonsetpos*20)-0.5,'c^',T(spwend*20),y1_6(spwend*20)-0.5,'go'); hold off; %
% xlabel('Time (s)'); ylabel('Voltage (mV)'); title('Raw LFP traces');

clear Fmax Fmin p A B

% figure;
% subplot(211);
% plot(t,x,'k'); xlim([0 t(end)]); % axis tight;
% hold on;
% plot(t(spwpeakpos),x(spwpeakpos),'ro', t(spwonsetpos),x(spwonsetpos),'c^'); hold off;% % 
% %ylabel('Voltage (mV)'); title('Filtered LFP traces ('num2str(length(spwpeakpos)) 'SWs (threshold = ' num2str(level) 'mV), 2-30 Hz');
% subplot(212);
% plot(t,y1-0.1,'k'); xlim([0 t(end)]); hold on;% axis tight; 
% plot(t(spwpeakpos),y1(spwpeakpos)-0.1,'ro', t(spwonsetpos),y1(spwonsetpos)-0.1,'c^'); hold off;
% xlabel('Time (s)'); ylabel('Voltage (mV)'); title('Raw LFP traces');

close(1000);
end

