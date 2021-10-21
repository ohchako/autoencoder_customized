%% Vm_fluctuation5.m %%
%% 3細胞記録以上に適用できるよう改変 %%
%% Update: 200729 %%
%膜電位が負および0にならないように検出するプログラム
% 脱分極・過分極を検出し、振幅の差が大きいほうをとる
% spwpeakposから-30～+40msのSWの時間index(spwallind2_2)

function [difamp1,allp,allf]=Vm_fluctuation5(L,fs,spwpeakpos,spwnum,t,y1,Vmall_2,yfilt,spwallind2_2);

% update 200729（上記プログラム＋過分極検出。振幅の差が大きいほうをとる）
for i=1:spwnum
    % 脱分極の検出
    [MAX(i,:),maxI(i,:)]=max(Vmall_2(i,:)); % maxの膜電位が一意に決まるのであれば１つのindex(maxI)が出る
    % maxの膜電位が２つ以上検出された場合、spwpeakに近い方をmaxの膜電位として検出する→knnsearchじゃなくて単純にspwpeakのindexと、検出されたindexの差分をとるのがよいかも（200729）
%     if size(maxI,2) >= 2
%         maxI(i,1)=knnsearch(spwallind2_2(I),spwallind2_2(i,601));
%     end
    f(i,1)=spwallind2_2(i,maxI(i,:)); % maxの膜電位がくるときのSPWのindex
    if  f(i,1)==spwallind2_2(i,1); % maxの膜電位が指定の範囲の左端にきた場合
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=epsp(i,1); % onsetとなる最小の膜電位はmaxの膜電位より20ms前
    elseif f(i,1)==spwallind2_2(i,end); % maxの膜電位が指定の範囲の右端にきた場合
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=epsp(i,end);% onsetとなる最小の膜電位はmaxの膜電位と同じ（実質ΔVmはゼロになる）
    else
        epsppeak(i,:)=[f(i)-round(0.2*L*fs):f(i)];
        epsp(i,:)=yfilt(epsppeak(i,:));
        mn(i,1)=min(epsp(i,:)); %maxの膜電位から-20msまでの範囲の中で最小の膜電位をとる
    end
    
    % 膜電位の変化を検出（脱分極）
    n(i,1)=min(find(epsp(i,:)==mn(i,1)));
    nn=epsppeak(i,:);
    nnn=n(i,:);
    ep(i,1)=(nn(nnn)); % index(プロット確認に使う)   
    epspamp(i,1)=MAX(i)-mn(i); %max-min＝膜電位変動 %peakから-20ms以内で見られる最小値をbaseとして膜電位変化を検出する
    
    % 過分極の検出
    [MIN(i,:),minI(i,:)]=min(Vmall_2(i,:));
   
   % minの膜電位が２つ以上検出された場合、spwpeakに近い方をmaxの膜電位として検出する
%     if size(minI,2) >= 2
%         minI(i,1)=knnsearch(spwallind2_2(minI),spwallind2_2(i,21));
%     end
    minf(i,1)=spwallind2_2(i,minI(i,:));
    
    % windowの左端にminが検出された場合、difampが0になるように
    if  minf(i,1)~=spwallind2_2(i,1); 
        ipspbase(i,:)=[minf(i)-round(0.2*L*fs):minf(i)];
        ipsp(i,:)=yfilt(ipspbase(i,:));
        mmn(i,1)=max(ipsp(i,:)); %min膜電位
    else
        ipspbase(i,:)=[minf(i)-round(0.2*L*fs):minf(i)];
        ipsp(i,:)=yfilt(ipspbase(i,:));
        mmn(i,1)=ipsp(i,1);
    end
    
    % 膜電位の変化を検出（過分極）
    minn(i,1)=min(find(ipsp(i,:)==mmn(i,1)));
    minnn=ipspbase(i,:);
    minnnn=minn(i,:);
    ip(i,1)=(minnn(minnnn)); % index(プロット確認に使う)
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

%膜電位差が正になるように。
% if any(difamp1==0);
%     ff=find(difamp1==0); %max-minがゼロになるものをfind。この行列を元のから削除し、新たなプログラムでつくった行列をhold onでプロットしてみる
%     l=4;
%     for i=1:numel(ff)
%         XX(i,:)=[spwallind2_2(ff(i,:),1):spwallind2_2(ff(i,:),end)];
%         VV(i,:)=[Vmall2_2(ff(i,:),:)];
%         XXX(i,:)=[XX(i,1):l:XX(i,end)];
%         zz(i,:)=spline(XX(i,:),VV(i,:),XXX(i,:));
%         Diff(i,:)=diff(zz(i,:)); 
%     end
%     
%     Diff=[Diff(:,1) Diff]; %Diffの1番目を最初の列に挿入し、1,2番目の列を同じ数にすることで差分がゼロになるようにする
%     Difflog=(Diff>0);
%     for i=1:numel(ff)
%         for j=1:size(Difflog,2)
%             if Difflog(i,j)==1
%                 Diffind(i,:)=j; %最初に1が入る各列のindexを返す。このあとに入るmax(peak)の値が膜電位変動のpeakの値。
%             %最終的にはmax(peak)から-20ms以内の最小値をmin(base)としたい！！
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
%         ccc(i,1)=max(y2sub((zzz(i,:):XX(i,end)))); %zzz=最小値のindex。cccはmax(peak)の値。
%         cc(i,1)=find(VV(i,:)==ccc(i,:)); %VVの行列の中で、何番目にmaxがくるか
%         MAX(i,1)=XX(i,cc(i,:)); %max(peak)のindex!
%         peakm20(i,:)=[MAX(i,:)-round(0.2*L*fs):1:MAX(i,:)]; %peakから-20ms以内で見られる最小値をbaseとして膜電位変化を検出する
%         Peak(i,:)=y2sub(peakm20(i,:));
%         Min(i,1)=min(Peak(i,:));%minの数値
%         minidx(i,1)=find(Peak(i,:)==Min(i,:)); %Peak行列の中で何番目にminがくるか
%         MIN(i,:)=peakm20(i,minidx(1,:));
%         difamp2(i,1)=ccc(i,1)-Min(i,1);
%     end
%     clear zzz ccc cc peakm20 Peak Min minidx
%     %次にdifamp2(difamp1ではゼロ）をdifamp1に挿入する。
%     difamp=[]; difamp=difamp1;
%     difamp(ff)=difamp2; %difamp=膜電位変化(mV)・・・求めたい行列！！
% 
%     f5=f; ep5=ep;
%     f5(ff,:)=[]; ep5(ff,:)=[];%max-minがゼロになるものを空行列に。
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
% 以下、ちゃんと検出できているかを確認するためのfigure
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



    
    
    