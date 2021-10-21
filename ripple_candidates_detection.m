%% ripple_candidates_detection.m
%% 18803 update, Nobuyoshi Matsumoto, OCU, 2018

%% USAGE
% EXAMPLE: ripple_detection_0('D:\data\nm0006\171208');
%
% Input:
% - basepath: atf�ۑ����Ă���t�H���_�B��F'D:\data\nm0006\171208'
%   basepath�ɐ}��v�Z���ʁi�����Excel�j���ۑ������B
%   basepath��atf�t�@�C����������Ȃ��ꍇ�͂����ǂݍ��ށB
%   basepath��atf�t�@�C������������ꍇ�́A���[�U�[�C���^�[�t�F�C�X����I�����ēǂݍ��ށB
% Output:
% datasaveoption��1�ɂ��Ă����ƁA�v�Z���ʁiAllCandidates�j�Ɛݒ�iAllSettings�j��basepath�ɕۑ������B
% figsaveoption��1�ɂ��Ă����ƁA�}�i.fig��.tif�j��basepath�ɕۑ������B
% xlssaveoption��1�ɂ��Ă����ƁA���b�v���^�m�C�Y�̐U�蕪���p��xls��basepath�ɕۑ������B

function ripple_candidates_detection(basepath)
tic;
cd(basepath);
clearvars -except basepath

%% Set patameters
ii=3; % import����atf�̃f�[�^�̉����LFP�̃f�[�^�������Ă��邩�Bdefault��3??

Fmax = 249; % �t�B���^���g���ő�l
Fmin = 120; % �t�B���^���g���ŏ��l
ThrSD = 3; % RMS��臒l�i�x�[�X���C�����牽SD���j)

rms_duration = 3; % in msec�@�@rms���v�Z���邽�߂�bin��duration
min_duration = 10; % in msec   ������������Ԃ��C�x���g�Ɣ��肵�Ȃ��B
min_isw_period = 5; % in msec ripple�Ԃ�������Z�����merge����B

datasaveoption = 1; % 1�̏ꍇ�A�ϐ��i�v�Z���ʁj��ۑ�����B
figsaveoption = 0; % 1�̏ꍇ�A�}�i�g���q .fig�E�E�E�d���I�j��ۑ�����B���̒i�K�ł͕ۑ����Ȃ��Ă��ǂ��i���b�v�����m��ł͂Ȃ����d���̂Łj�B
tifsaveoption = 1; % 1�̏ꍇ�Atiff�̐}��ۑ�����B.fig�t�@�C�����͌y���B
xlsmakeoption = 1; % 1�̏ꍇ�AExcel�t�@�C�����쐬����B

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
lfp0=d(:,ii); %LFP���܂܂�����w��iLFP�̔g�`�j
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
y=filtfilt(B,A,lfp0); % filter���LFP�̔g�`�i���ƂŁAy��lfpf�Ɩ��O��ς���j

%% Detect ripple candidates (onset & offset): calculate RMS, threshold, zscore
rmsdurainidx=round(rms_duration*fs/1000); % rms_duration�̓C���f�b�N�X�Ɋ��Z����Ƃ������ɂȂ邩�B�i3��4�Ȃǂ̐����ЂƂ̃X�J���[���o�Ă���j
forenlarge_y=zeros(rmsdurainidx,1);
forenlarge_y(:)=mean(y);
temp_y=[forenlarge_y;y;forenlarge_y];
temp_y=temp_y-mean(temp_y); % baseline��0�ɂ����B
rms_y=zeros(size(y));

for i=1:size(y(:,1))   %��敽�ϕ������CRMS
   rms_y(i)=sqrt(sum(temp_y(i:i+rmsdurainidx-1).^2)/rmsdurainidx);   
end

mean_rms_y=zeros(size(rms_y)); mean_rms_y(:)=mean(rms_y);
std_rms_y=zeros(size(rms_y)); std_rms_y(:)=std(rms_y);

zrippledetection=(rms_y-mean_rms_y)./std_rms_y;

rippleidx=zeros(size(lfp0)); % ���b�v���̌��͂�����1�����Ă����i�ŏ���zeros�őS��0�ɂ��Ă����j
% lfp0����x�N�g���Ȃ�Arippleidx����x�N�g���ɂȂ�B

rippleidx(zrippledetection>=ThrSD)=1; % threshold�𒴂���������1������B

%% Remove ripple candidates if each duration is very short
% ���b�v����duration���Z��������̂�0�ɂ���imin_duration���g���j
diffrippleidx=[0;diff(rippleidx)]; % �����Ƃ��������B��x�N�g���ɂȂ�͂��B
onsetidx =  find(diffrippleidx== 1); % onset�̌��
offsetidx= (find(diffrippleidx==-1))-1; % offset�̌��

numon = numel(onsetidx);
numoff= numel(offsetidx);
% ��numon��numoff�͓������ł����āA���Aonsetidx�͏��offsetidx�����������Ȃ���΂Ȃ�Ȃ��̂ŁA
% �ȉ��̂悤�ɂ���B
if numon > numoff
    onsetidx(end)=[];
elseif numon < numoff
    offsetidx(1)=[];
end

rippleinfo = [[1:numel(onsetidx)]', onsetidx, offsetidx,...
              onsetidx/fs, offsetidx/fs];
% 1���ripple�ԍ��A2���ripple onset��index�A3���ripple offset��index�A
% 4���ripple onset�̎����i�b�j�A5���ripple offset�̎����i�b�j�B

for i=1:size(rippleinfo,1) % �Z�����郊�b�v���͍폜�i�܂�0�Œu�������āAfor�����I�������폜����j
    if (rippleinfo(i,5)-rippleinfo(i,4)) <= min_duration/1000 % msec��sec�ɂ��邱�Ƃɒ���
        rippleinfo(i,:)=[0 0 0 0 0];
    end
end
i=rippleinfo(:,1); excludeidx=find(i==0);
rippleinfo(excludeidx,:)=[];
rippleinfo(:,1)=[1:size(rippleinfo,1)]'; % ripple�ԍ���U�蒼�����i1,2,3,4...�ɂȂ�悤�ɂ����j
rippleinfo=[rippleinfo,zeros(size(rippleinfo,1),3)]; % 6���ripple��duration,7���peak��index, 8���peak�̎���������
rippleinfo(:,6)=(rippleinfo(:,5)-rippleinfo(:,4))*1000; % duration�A�P�ʂ́u�~���b�v

%% Merge ripple candidates if inter-ripple intarvals are very short
% �߂����郊�b�v����merge����
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
for i=1:size(rippleinfo,1) % �s�[�N�̌��o
    eachrms=rms_y(rippleinfo(i,2):rippleinfo(i,3));
    pkidxcand=find(eachrms==max(eachrms),1,'first');
    rippleinfo(i,7)=rippleinfo(i,2)+pkidxcand-1; %-1�����Ȃ��ƃs�[�N�̃C���f�b�N�X����ЂƂ���Ă��܂�
    rippleinfo(i,8)=time(rippleinfo(i,7));
end

rippleinfo(:,4)=time(rippleinfo(:,2));
rippleinfo(:,5)=time(rippleinfo(:,3));

%% Draw figures (and save if required)
% �`�� 
SS = get(0, 'ScreenSize'); % �X�N���[���T�C�Y���擾
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
% �����ŁAripple�̐U�蕪���i�m�C�Y or �V�O�i���j�p��xls�t�@�C�������B
if xlsmakeoption
    xlsinfo=[rippleinfo(:,1), zeros(size(rippleinfo,1),1),...
             rippleinfo(:,4), rippleinfo(:,8), rippleinfo(:,5)];
         % 1���ripple�ԍ��A2���0 or 1 or 2����������ł����B
         % 3���ripple onset�̎����i�b�j�A4���ripple peak�̎����i�b�j�A
         % 5���ripple offset�̎����i�b�j�B
    xlswrite(xlsname,xlsinfo); % xlsx�t�@�C���������ō��B
    % ����xls�t�@�C����2��ɁA���l�����Ă����B
    % 0: �m�C�Y
    % 1: �P�ƂŃ��b�v��
    % 2: �Ȃ��郊�b�v���̍ŏ�
    % 0: �Ȃ��郊�b�v���̍ŏ��ł��Ō�ł��Ȃ��Ƃ���i����͂����Ă��Ȃ��Ă������j
    % 3: �Ȃ��郊�b�v���̍Ō�i2�Ɠ������ɂȂ�͂��j
end

%% Save candidates
% �ϐ���ݒ��ۑ�����
cd(basepath);
lfpf = y; % �o���₷�����邽�߂ɁAfilter���lfp�̕ϐ��̖��O��lfpf�Ƃ��������Bfiltered LFP�Ƃ����Ӗ���lfpf�Ƃ����B
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
% �^�C�g������Text�I�u�W�F�N�g�ɃA���_�[�o�[�u_�v��n�b�g�u^�v��\���������̂ł����A
% ���ꂼ��̋L���̒���̕��������t?�������Ə�t�������ɂȂ��Ă��܂��܂��B
% https://jp.mathworks.com/matlabcentral/answers/98127-text-_
%
% find��first��last�ƕ��p����
% https://jp.mathworks.com/help/matlab/ref/find.html
% k = find(X,n,direction) �́Adirection �� 'last' �̏ꍇ�A
% X ���̔�[���v�f�ɑΉ�����C���f�b�N�X���Ōォ�琔���� n ���o���܂��B
% direction �̊���l�� 'first' �ł���A��[���v�f�ɑΉ�����C���f�b�N�X���ŏ����琔���� n ���o���܂��B
