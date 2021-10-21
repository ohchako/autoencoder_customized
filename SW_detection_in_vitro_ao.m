%% SW_detection_in_vitro_nm5,m %%
%% 20160817 filter_spw_peak_num.m������
%% to find local maxima of x (>level). ��ɓ˂��o����peak��position��amplitude, ������������
%% �s�[�N��臒l�ƃx�[�X���C���̕␳��mfile��ōs��
%% data reduction���sampling frequency��1000 Hz�𐄏�
%% y1: local field potentials (mV)
%% Update: 20170522 %%

function [spwpeakpos, spwpeaktime, spwpeakamp, spwnum, spwfrq, spwonsetpos, spwdur, spwend, y1sw]=SW_detection_in_vitro_ao(y1,t,T,y1_6);
%[spwpeakpos, spwpeaktime, spwpeakamp, spwnum, spwfrq, spwonsetpos, spwdur, spwend, y1sw, x]=ao_SW_detection_in_vitro(y1,t,T,y1_6);
% spwpeakamp (��V)
% spwpeaktime (s)
% spwfrq (Hz)

%%%%%%%%%%%%%
%�p�����^�ݒ�

p=0;%baseline��␳����
fs=1/(t(2)-t(1)); % sampling frequency (Hz)
Fmax=30; %�t�B���^���g���ő�l
Fmin=2; %�t�B���^���g���ŏ��l
%%%%%%%%%%%%%
%filter�������܂�
%%%%%%%%%%%%%%%%%%%%%%%%%%

[B,A]=butter(3,[Fmin/(fs*0.5) Fmax/(fs*0.5)]);
x=filtfilt(B,A,y1); y1sw=x;

% figure;
% subplot(211);
% plot(t,y1,'k');
% subplot(212);
% plot(t,x,'k');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%�S�̂�SD�l�A�x�[�X���C���m�C�Y��SD�l���Z�o���܂��B
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
stdv=std(x) %SD
stdv=2*std(x) %2SD
stdv=3*std(x) %3SD
stdv=4*std(x) %4SD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%baseline�����߂�
xmin=25;%plot����͈͂�K���ɂ��߂�
xmax=35;
ymin=-0.01;
ymax=0.1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1000); plot(t,x,'k');
axis([xmin, xmax, ymin, ymax]);
baseline=ginput(2);
basexmin=baseline(1,1);
basexmax=baseline(2,1);

basewave=x([round(basexmin*fs):round(basexmax*fs)],:);%�����l�łȂ��Ă͂Ȃ�܂���
SD=std(basewave);
ave=mean(basewave);
ave+SD
ave+2*SD
ave+3*SD
ave+4*SD
ave+5*SD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
level=input('level(mV) >>\n');%臒l�̐ݒ�i��V�j

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
% spwpeakpos=idx(pidx); %�C���f�b�N�X�����߂�
spwpeakpos=spwpeakpos(spwpeakpos>0.1*fs); spwpeakpos=spwpeakpos(spwpeakpos<(round(t(end))-0.1)*fs);
% ��nobumossy�̒��ŁApeak��0.1�b�O����0.1�b��܂ł�trace���d�˂邩��A�ŏ��ƍŌ�̕��̃s�[�N�͕s�v�B
spwpeakamp=x(spwpeakpos)*1000; % spwpeakamp�̓�V�P�ʂɂ���
% spwpeaktime=spwpeakpos/fs;
spwnum=size(spwpeakpos,1);
spwfrq=size(spwpeakpos,1)/round(t(end));
spwpeaktime=t(spwpeakpos);


spwonsetpos=zeros(size(spwpeakpos)); %SPW onset��index

for j=1:numel(spwpeakpos); % �s�[�N����50ms�O�܂ł�filtered trace�ɂ����āAbaseline+sd�𒴂����_��onset�Ƃ���
    AA=[round(spwpeakpos(j)-0.1*fs):round(spwpeakpos(j))]; %�s�[�N����50ms�O�܂ł̃C���f�b�N�X�i50ms�ɂ����͓̂K�������A���͈͓̔���onset��������ꂻ���j
    B=find(y1sw(AA)>(ave+SD));
    if size(B,2)>1 % �s�x�N�g�����x�N�g���ɂ���
        B=B';
    end
    dB=[2;diff(B)]; % �ŏ��̐�����1���傫�����ɂ��Ă����΁Aif isequal(dB,ones(size(dB))==1�Ƃ������Ȃ��Ă悢
    c=max(find(dB>1)); % �s�[�N�Ɉ�ԋ߂�������肽���̂ŁAmax����ꂽ
    BB=B(c);
    CC=AA(BB); % ���ꂪfiltered trace�ɂ�����SD�𒴂������߂Ă̎�����index
    
    % �s�[�N�̒��O�̃g���t���܂߂�ave+sd�𒴂��Ă��܂����ꍇ�͎d���Ȃ��̂ŁA�s�[�N�̒��O�̃g���t����SD�����������_���Ƃ��Ă���
    yo=-y1sw(CC:AA(end)); % �g�`�𔽓]�����ăg���t��������
    [~,postr]=findpeaks(yo,'MINPEAKDISTANCE',5);% round(0.005*fs));��5ms�iround�ł͂��܂������Ȃ��̂ŁA���l�𒼐ړ����j % trough�F�g���t % �g���t�̑O��10ms�ŗǂ����������Ƃ肠����5ms�ɂ��Ă݂�
    if isempty(postr)==0 % ���s�[�N���O�̃g���t���܂�ave+sd�𒴂����ꍇ
        cc=max(postr);
        AAasc=AA(cc:end); % cc����͑�������悤�ɂȂ��Ă���
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
for j=1:numel(spwpeakpos); % �s�[�N����50ms��܂ł�filtered trace�ɂ����āApeakpos����baseline+SD�ɖ߂����_��end
    A2=[round(spwpeakpos(j)):round(spwpeakpos(j)+0.1*fs)]; %�s�[�N����100ms��܂ł̃C���f�b�N�X�i100ms�ɂ����͓̂K�������A���͈͓̔���end��������ꂻ���j
    B2=find(y1sw(A2)<=(ave+SD));
    c2=min(B2);% �s�[�N�Ɉ�ԋ߂�������肽���̂ŁAmin����ꂽ�B�܂�s�[�N�ɍł��߂�ave+SD�������index�B%�����܂ł���
    BB2=A2(c2); %ave�������index�̒l
%     CC2=A2(BB2); % ���ꂪfiltered trace�ɂ�����SD�𒴂������߂Ă̎�����index
%   spwend(j,1)=BB2;
end

spwdur=zeros(size(spwpeakpos));
for k=1:numel(spwpeakpos)
spwdur(k,1)=spwend(k,1)-spwonsetpos(k,1);
end
spwdur=spwdur/fs; %�����is�j
clear k

%�ȉ���figure���m�F�p�ɕ`�}����ۂ�17�|26�s�����s���A�K�v�ȕϐ�(x)�������B
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

