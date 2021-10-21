%% Vm_fluctuation5.m %%
%% 3�זE�L�^�ȏ�ɓK�p�ł���悤���� %%
%% Update: 200729 %%
%���d�ʂ��������0�ɂȂ�Ȃ��悤�Ɍ��o����v���O����
% �E���ɁE�ߕ��ɂ����o���A�U���̍����傫���ق����Ƃ�
% spwpeakpos����-30�`+40ms��SW�̎���index(spwallind2_2)

function [difamp1,allp,allf]=Vm_fluctuation5(L,fs,spwpeakpos,spwnum,t,y1,Vmall_2,yfilt,spwallind2_2);

% update 200729�i��L�v���O�����{�ߕ��Ɍ��o�B�U���̍����傫���ق����Ƃ�j
for i=1:spwnum
    % �E���ɂ̌��o
    [MAX(i,:),maxI(i,:)]=max(Vmall_2(i,:)); % max�̖��d�ʂ���ӂɌ��܂�̂ł���΂P��index(maxI)���o��
    % max�̖��d�ʂ��Q�ȏ㌟�o���ꂽ�ꍇ�Aspwpeak�ɋ߂�����max�̖��d�ʂƂ��Č��o���遨knnsearch����Ȃ��ĒP����spwpeak��index�ƁA���o���ꂽindex�̍������Ƃ�̂��悢�����i200729�j
%     if size(maxI,2) >= 2
%         maxI(i,1)=knnsearch(spwallind2_2(I),spwallind2_2(i,601));
%     end
    f(i,1)=spwallind2_2(i,maxI(i,:)); % max�̖��d�ʂ�����Ƃ���SPW��index
    if  f(i,1)==spwallind2_2(i,1); % max�̖��d�ʂ��w��͈̔͂̍��[�ɂ����ꍇ
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=epsp(i,1); % onset�ƂȂ�ŏ��̖��d�ʂ�max�̖��d�ʂ��20ms�O
    elseif f(i,1)==spwallind2_2(i,end); % max�̖��d�ʂ��w��͈̔͂̉E�[�ɂ����ꍇ
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=epsp(i,end);% onset�ƂȂ�ŏ��̖��d�ʂ�max�̖��d�ʂƓ����i������Vm�̓[���ɂȂ�j
    else
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=min(epsp(i,:)); %max�̖��d�ʂ���-20ms�܂ł͈̔͂̒��ōŏ��̖��d�ʂ��Ƃ�
    end
    
    % ���d�ʂ̕ω������o�i�E���Ɂj
    n(i,1)=min(find(epsp(i,:)==mn(i,1)));
    nn=epsppeak(i,:);
    nnn=n(i,:);
    ep(i,1)=(nn(nnn)); % index(�v���b�g�m�F�Ɏg��)   
    epspamp(i,1)=MAX(i)-mn(i); %max-min�����d�ʕϓ� %peak����-20ms�ȓ��Ō�����ŏ��l��base�Ƃ��Ė��d�ʕω������o����
    
    % �ߕ��ɂ̌��o
    [MIN(i,:),minI(i,:)]=min(Vmall_2(i,:));
   
   % min�̖��d�ʂ��Q�ȏ㌟�o���ꂽ�ꍇ�Aspwpeak�ɋ߂�����max�̖��d�ʂƂ��Č��o����
%     if size(minI,2) >= 2
%         minI(i,1)=knnsearch(spwallind2_2(minI),spwallind2_2(i,21));
%     end
    minf(i,1)=spwallind2_2(i,minI(i,:));
    
    % window�̍��[��min�����o���ꂽ�ꍇ�Adifamp��0�ɂȂ�悤��
    if  minf(i,1)~=spwallind2_2(i,1); 
        ipspbase(i,:)=[minf(i)-round(0.2*L*fs):minf(i)];
        ipsp(i,:)=yfilt(ipspbase(i,:));
        mmn(i,1)=max(ipsp(i,:)); %min���d��
    else
        ipspbase(i,:)=[minf(i)-round(0.2*L*fs):minf(i)];
        ipsp(i,:)=yfilt(ipspbase(i,:));
        mmn(i,1)=ipsp(i,1);
    end
    
    % ���d�ʂ̕ω������o�i�ߕ��Ɂj
    minn(i,1)=min(find(ipsp(i,:)==mmn(i,1)));
    minnn=ipspbase(i,:);
    minnnn=minn(i,:);
    ip(i,1)=(minnn(minnnn)); % index(�v���b�g�m�F�Ɏg��)
    ipspamp(i,1)=MIN(i)-mmn(i);
    
    if abs(epspamp(i,1)) > abs(ipspamp(i,1))
        difamp1(i,1)=epspamp(i,1);
        allp(i,1)=ep(i,1);
        allf(i,1)=f(i,1);
    else
        difamp1(i,1)=ipspamp(i,1);
        allp(i,1)=ip(i,1);
        allf(i,1)=minf(i,1);
    end
    
end

 
% figure; plot(t,y2sub,'k',t(f),y2sub(f),'ro',t(ep),y2sub(ep),'ys','MarkerSize',10); 
% hold on; plot(t,40*y1-40);

%���d�ʍ������ɂȂ�悤�ɁB
% if any(difamp1==0);
%     ff=find(difamp1==0); %max-min���[���ɂȂ���̂�find�B���̍s������̂���폜���A�V���ȃv���O�����ł������s���hold on�Ńv���b�g���Ă݂�
%     l=4;
%     for i=1:numel(ff)
%         XX(i,:)=[spwallind2_2(ff(i,:),1):spwallind2_2(ff(i,:),end)];
%         VV(i,:)=[Vmall2_2(ff(i,:),:)];
%         XXX(i,:)=[XX(i,1):l:XX(i,end)];
%         zz(i,:)=spline(XX(i,:),VV(i,:),XXX(i,:));
%         Diff(i,:)=diff(zz(i,:)); 
%     end
%     
%     Diff=[Diff(:,1) Diff]; %Diff��1�Ԗڂ��ŏ��̗�ɑ}�����A1,2�Ԗڂ̗�𓯂����ɂ��邱�Ƃō������[���ɂȂ�悤�ɂ���
%     Difflog=(Diff>0);
%     for i=1:numel(ff)
%         for j=1:size(Difflog,2)
%             if Difflog(i,j)==1
%                 Diffind(i,:)=j; %�ŏ���1������e���index��Ԃ��B���̂��Ƃɓ���max(peak)�̒l�����d�ʕϓ���peak�̒l�B
%             %�ŏI�I�ɂ�max(peak)����-20ms�ȓ��̍ŏ��l��min(base)�Ƃ������I�I
%             break
%             end
%         end
%     end
%     
% %     D=find(Diffind==0);
% %     Diffind(D)=[];
%     
%     clear zz Diff Difflog
%     % XXXnumel=size(XXX,2); VVnumel=size(VV,2);
%     % XVnumel=round(VVnumel/XXXnumel);
%     for i=1:numel(ff)
%         zzz(i,:)=XXX(i,Diffind(i));
%         ccc(i,1)=max(y2sub((zzz(i,:):XX(i,end)))); %zzz=�ŏ��l��index�Bccc��max(peak)�̒l�B
%         cc(i,1)=find(VV(i,:)==ccc(i,:)); %VV�̍s��̒��ŁA���Ԗڂ�max�����邩
%         MAX(i,1)=XX(i,cc(i,:)); %max(peak)��index!
%         peakm20(i,:)=[MAX(i,:)-round(0.2*L*fs):1:MAX(i,:)]; %peak����-20ms�ȓ��Ō�����ŏ��l��base�Ƃ��Ė��d�ʕω������o����
%         Peak(i,:)=y2sub(peakm20(i,:));
%         Min(i,1)=min(Peak(i,:));%min�̐��l
%         minidx(i,1)=find(Peak(i,:)==Min(i,:)); %Peak�s��̒��ŉ��Ԗڂ�min�����邩
%         MIN(i,:)=peakm20(i,minidx(1,:));
%         difamp2(i,1)=ccc(i,1)-Min(i,1);
%     end
%     clear zzz ccc cc peakm20 Peak Min minidx
%     %����difamp2(difamp1�ł̓[���j��difamp1�ɑ}������B
%     difamp=[]; difamp=difamp1;
%     difamp(ff)=difamp2; %difamp=���d�ʕω�(mV)�E�E�E���߂����s��I�I
% 
%     f5=f; ep5=ep;
%     f5(ff,:)=[]; ep5(ff,:)=[];%max-min���[���ɂȂ���̂���s��ɁB
% 
%     figure(1);
%     plot(t,y2sub,'k',...
%         t(f5),y2sub(f5),'ro',...
%         t(ep5),y2sub(ep5),'ys','MarkerSize',10);
%     hold on; plot(t,40*y1-40);
%     hold on; plot(t(MIN),y2sub(MIN),'co',t(MAX),y2sub(MAX),'mo');
% 
% else
% 
% �ȉ��A�����ƌ��o�ł��Ă��邩���m�F���邽�߂�figure
% A=yfilt(:,1); AA=yfilt(:,2); AAA=yfilt(:,3);
% B=f(:,1); BB=f(:,2); BBB=f(:,3);
% C=ep(:,1); CC=ep(:,2); CCC=ep(:,3);
%     figure;
%     plot(t,A+40,'k',...
%         t(B),A(B)+40,'ro',...
%         t(C),A(C)+40,'gs','MarkerSize',5); hold on;
%     plot(t,AA+20,'k',...
%         t(BB),AA(BB)+20,'ro',...
%         t(CC),AA(CC)+20,'gs','MarkerSize',5); hold on;
%     plot(t,AAA,'k',...
%         t(BBB),AAA(BBB),'ro',...
%         t(CCC),AAA(CCC),'gs','MarkerSize',5); hold on;
%     plot(t,70*y1); hold off;      
    

end



    
    
    