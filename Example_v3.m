% Ref. Antonio Cicone, Haomin Zhou. "Multidimensional Iterative Filtering method
%      for the decomposition of high-dimensional non-stationary signals".
%      Cambridge Core in Numerical Mathematics: Theory, Methods and
%      Applications, Volume 10, Issue 2, Pages 278-298, 2017.
%      doi:10.4208/nmtma.2017.s05
%
%      Stefano Sfarra, Antonio Cicone, Bardia Yousefi, Stefano Perilli,
%      Leonardo Robol, Xavier P.V. Maldague.
%      "Maximizing the detection of thermal imprints in civil engineering
%      composites after a thermal stimulus - The contribution of an
%      innovative mathematical pre-processing tool: the 2D Fast Iterative
%      Filtering algorithm. Philosophy, comparisons, numerical, qualitative
%      and quantitative results". 2021. Submitted
%

%%

z2=[];

iMax=8;
for i=0:iMax
    T=2-2*i/10;
    t=0:0.02:T;
    z21=t/T;
    z22=fliplr(z21);
    if not(i==iMax) && not(i==0)
        z2t=[z21 z22(2:end-1)];
    elseif i==iMax
        z2t=[z21 z22(2:end)];
    else
        z2t=[z21(ceil(end/2):end) z22(2:end-1)];
    end
    z2=[z2 z2t];
end
z2t=fliplr(z2);
z2=[z2 z2t(2:end)];
z2((length(z2)+1)/2)=0;
z2=z2-0.5;

%%
f1=plot_2D_v1(z2((length(z2)+1)/2:end),200);
f1=f1/max(max(f1))*2.5;
surf(f1,'edgecolor','none')
%%

t2=linspace(0,10,length(z2));
z3=sin(pi*0.7*t2);

%%

f2=plot_2D_v1(z3((length(z3)+1)/2:end),200);
f2=f2/max(max(f2))*5;
surf(f2,'edgecolor','none')

%%

x=linspace(0,10+1/2,200*2+1);
f3=0.5*x.'*ones(1,length(x));
surf(f3,'edgecolor','none')

%%

f=f1+f2+f3;
%%

figure
h=surf(f);
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(f,1) 1 size(f,2) floor(min(min(f))) ceil(max(max(f)))])

%% Section of the signal in the middle

figure
plot((f(:,(end+1)/2)),'k','Linewidth',2)
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(f,1) floor(min(min(f))) ceil(max(max(f)))])

%% IMFs
close all
opts=Settings_FIF2_v2('delta',0.01,'NIMFs',1,'alpha','ave');
tic
[IMF_temp1,SDlog1] = FIF2_v3(f,opts);
toc
opts=Settings_FIF2_v2('delta',0.1,'NIMFs',1,'alpha',90);
tic
[IMF_temp2,SDlog2] = FIF2_v3(IMF_temp1(:,:,end),opts);
toc

IMF(:,:,1)=IMF_temp1(:,:,1);
IMF(:,:,2:3)=IMF_temp2(:,:,1:2);

%%

figure
h=surf(IMF(:,:,1));
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(IMF(:,:,1),1) 1 size(IMF(:,:,1),2) floor(min(min(IMF(:,:,1)))) ceil(max(max(IMF(:,:,1))))])


%% Section of IMF1 in the middle

figure
plot((f1(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF(:,(end+1)/2,1)),'k','Linewidth',2)
legend('Ground truth','IMF_1')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF(:,:,1),1) floor(min(min(IMF(:,:,1)))) ceil(max(max(IMF(:,:,1))))])


%% Second IMF

fprintf('\n\n IMF 2\n\n')


figure
h=surf(IMF(:,:,2));
set(h, 'edgecolor','none')
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
colorbar
set(gca,'fontsize', 25);
axis([1 size(IMF(:,:,2),1) 1 size(IMF(:,:,2),2) floor(min(min(IMF(:,:,2)))) ceil(max(max(IMF(:,:,2))))])


%% Section of IMF2 in the middle

figure
plot((f2(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF(:,(end+1)/2,2)),'k','Linewidth',2)
legend('Ground truth','IMF_2')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF(:,:,2),1) floor(min(min(IMF(:,:,2)))) ceil(max(max(IMF(:,:,2))))])

%% Section of the Remainder along the diagonal

figure
plot((f3(:,(end+1)/2)),'r--','Linewidth',2)
hold on
plot((IMF(:,(end+1)/2,3)),'k','Linewidth',2)
legend('Ground truth','Remainder')
set(gca,'fontsize', 25);
set(gcf,'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
axis([1 size(IMF(:,:,3),1) floor(min(min(IMF(:,:,3)))) ceil(max(max(IMF(:,:,3))))])



%% Difference between IMF1 and the GT first component

h=surf(IMF(:,:,1)-f1);
set(h, 'edgecolor','none')


%% Difference between IMF2 and the GT second component

h=surf(IMF(:,:,2)-f2);
set(h, 'edgecolor','none')







