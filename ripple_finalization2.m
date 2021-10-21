%% ripple_finalization2.m
%% 190108 update, Nobuyoshi Matsumoto, Univ Tokyo, 2019
% ripple_finalization.m����onset��offset�������Ȃ̂�2SD�Ŏ�邱�Ƃɂ����B
% ripple_finalization.m�̌�Ɏg���B
% �ۑ������t�@�C���̂����A�t�@�C������final2�ƂȂ��Ă�����̂��A���̃v���O�����ŋ��߂��ϐ���}

%% USAGE
% EXAMPLE: ripple_finalization2('D:\data\nm0006\171208');
%
% Input:
% - basepath: ripple or noise��ڂŐU�蕪����Excel��ۑ����Ă���t�H���_
%   basepath�ɐ}��v�Z���ʂ��ۑ������B
%   basepath��mat�t�@�C����������Ȃ��ꍇ�͂����ǂݍ��ށB
%   basepath��atf�t�@�C������������ꍇ�́A���[�U�[�C���^�[�t�F�C�X����I�����ēǂݍ��ށB
% Output:
% datasaveoption��1�ɂ��Ă����ƁA�v�Z���ʁiAllResults�j�Ɛݒ�iAllSettings�j��basepath�ɕۑ������B
% figsaveoption��1�ɂ��Ă����ƁA�}�i.fig��.tif�j��basepath�ɕۑ������B
% deleteoldresults��1�ɂ��Ă����ƁAripple_candidates_detection.m�Ōv�Z���ꂽ�Â����ʂ���������B

function ripple_finalization2(basepath)
tic; cd(basepath); clearvars -except basepath

%% Set patameters
datasaveoption = 1; % 1�̏ꍇ�A�ϐ��i�v�Z���ʁj��ۑ�����B
figsaveoption = 0; % 1�̏ꍇ�A�}�i.fig�j��ۑ�����B����̏ꍇ�A.tif���d���I
tifsaveoption = 1; % 1�̏ꍇ�Atiff�̐}��ۑ�����B.fig�t�@�C�����͌y���B
deleteoldresults = 0; % 1�̏ꍇ�Aripple_candidates_detection.m�Ōv�Z���ꂽ�Â����ʂ���������B

ThrONOFF = 2; % �I���Z�b�g�ƃI�t�Z�b�g�����߂邽�߂�臒l
FminSGamma = 20; % �X���[�K���}�ш�̉������g��
FmaxSGamma = 40; % �X���[�K���}�ш�̏�����g��

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
ARlfp0=d(:,AllSettings.ii); %LFP���܂܂�����w��iLFP�̔g�`�j
fs=ARfs; % sampling frequency
clear d

% -- Filter LFP traces
[B,A]=butter(3,[AllSettings.Fmin/(fs*0.5) AllSettings.Fmax/(fs*0.5)]);
ARlfpf=filtfilt(B,A,ARlfp0);

% -- Calculate RMS
rmsdurainidx=round(AllSettings.rms_duration*fs/1000);
forenlarge_y=zeros(rmsdurainidx,1);
forenlarge_y(:)=mean(ARlfpf); % ���g�`�Ƀt�B���^�[���������g�`�̓d�ʂ̕��ϒl���c��rmsdurationidx���ׂ��B
temp_y=[forenlarge_y;ARlfpf;forenlarge_y]; % ���[��␳���邽�߂ɁA�ŏ��ƍŌ��forenlarge_y�����������B
temp_y=temp_y-mean(temp_y); % baseline��0�ɂ����Btemp_y����mean�������Ă��A���temp_y�𕽋ς��āA�O���0�̗���������Ă����ʂ͓����B
ARrms_y=zeros(size(ARlfpf));

for i=1:size(ARlfpf(:,1))   %��敽�ϕ������CRMS
   ARrms_y(i)=sqrt(sum(temp_y(i:i+rmsdurainidx-1).^2)/rmsdurainidx);   
end
clear i

Zripplerms=(ARrms_y-mean(ARrms_y)*ones(size(ARrms_y)))./(std(ARrms_y)*ones(size(ARrms_y)));

%% Detect ripple peaks: find max RMS
% peak�����ߒ����K�v�͂Ȃ��B

%% Detect ripple onsets and offsets
ONOFFidx=zeros(size(ARlfp0)); % ���b�v���̌��͂�����1�����Ă����i�ŏ���zeros�őS��0�ɂ��Ă����j
% ARlfp0����x�N�g���Ȃ�Arippleidx����x�N�g���ɂȂ�B
ONOFFidx(Zripplerms>=ThrONOFF)=1; % threshold�𒴂���������1������B
DiffONOFFidx=[ONOFFidx(1);diff(ONOFFidx)];
ripplefinalindex2=AllResults.ripplefinalindex; % onset��offset�����ߒ����A����5��ɁA�s�[�N����SD����������
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
             Zripplerms(ripplefinalindex2(:,3))]; % ���̎��_�ł�14��ɂ��Ă����B
%  1��@���b�v���ԍ�
%  2��@�I���Z�b�g�̃C���f�b�N�X
%  3��@�s�[�N�̃C���f�b�N�X
%  4��@�I�t�Z�b�g�̃C���f�b�N�X
%  5��@���b�v���ԍ�
%  6��@�I���Z�b�g�̎����i�b�j
%  7��@�s�[�N�̎����i�b�j
%  8��@�I�t�Z�b�g�̎����i�b�j
%  9��@�������ԁi�b�j
% 10��@Symmetry index (=(a-b)/(a+b): a �s�[�N�ƃI���Z�b�g�̋����Ab �I�t�Z�b�g�ƃs�[�N�̋����j
% 11��@���b�v���ш�̃E�F�[�u���b�g�p���[�i�E�F�[�u���b�g�X�y�N�g���̋��x�̘a�����ԂŋK�i���AFmin�`Fmax Hz�j
%   �@�@�@�c����Ă��鎞�ԑт��Z������FFT���ƃX�y�N�g�����e������̂ŁA�E�F�[�u���b�g�ɂ���
% 12��@�X���[�K���}�ш��FFT�p���[�i�E�F�[�u���b�g�X�y�N�g���̋��x�̘a�����ԂŋK�i���A20�`40 Hz�j
% 13��@���b�v���ш�̃s�[�N���g���i�E�F�[�u���b�g�X�y�N�g���̃��b�v���ш�ɂ����āAPSD���ő�ɂȂ���g���j
% 14��@���b�v���̃s�[�N����SD��

f_axis=[1:.01:2.5];
r1 =find(f_axis>log10(AllSettings.Fmin),1,'first'); % ���b�v���̍Œ���g���̏�p�ΐ���f_axis�ɂ�����u�C���f�b�N�X�v
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
% �`�� 
SS = get(0, 'ScreenSize'); % �X�N���[���T�C�Y���擾
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
% �A���_�[�o�[�����̂܂܃^�C�g���ɕ\���ł���悤�ɂȂ�
%�i��������Ȃ��ƁA�A���_�[�o�[�̎��̕��������t���ɂȂ�j
% �Ō��Reference���Q�l�ɁB
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
if deleteoldresults % �ŏ��̐ݒ�ŁAdeleteoldresults��1�ɂ��Ă����ƁA�Â��ϐ���}�������ł���B
    delete([ARdataname(1:end-4),'.mat']);
    delete([ARdataname(1:end-4),'.tif']);
    delete([ARdataname(1:end-4),'.fig']);
end

%% Save results and settings
% �ϐ���ݒ��ۑ�����
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
% �^�C�g������Text�I�u�W�F�N�g�ɃA���_�[�o�[�u_�v��n�b�g�u^�v��\���������̂ł����A
% ���ꂼ��̋L���̒���̕��������t?�������Ə�t�������ɂȂ��Ă��܂��܂��B
% https://jp.mathworks.com/matlabcentral/answers/98127-text-_
%
% find��first��last�ƕ��p����
% https://jp.mathworks.com/help/matlab/ref/find.html
% k = find(X,n,direction) �́Adirection �� 'last' �̏ꍇ�A
% X ���̔�[���v�f�ɑΉ�����C���f�b�N�X���Ōォ�琔���� n ���o���܂��B
% direction �̊���l�� 'first' �ł���A��[���v�f�ɑΉ�����C���f�b�N�X���ŏ����琔���� n ���o���܂��B
