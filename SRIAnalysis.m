%
%
%
% Setin up the basin info
load('InputDataRunoff.mat') % Runoff from the models
BasinInd=[168,358,195,124,58,353,16,63,54,393,71,352,...
88,351,177,354,107,392,259,336,69,19,81,47,41,...
51,126,82,25,49,187,118,83,128,80,110,59,149,6,...
23,344,292,123,131,399,67,346,199,345,307,350,270,...
24,235,65,10];
BasinNam={'INDUS','TARIM','BRAHMAPUTRA','ARAL SEA','COPPER','GANGES','YUKON','ALSEK','SUSITNA','BALKHASH','STIKINE','SANTA CRUZ',...
'FRASER','BAKER','YANGTZE','SALWEEN','COLUMBIA','ISSYK-KUL','AMAZON','COLORADO','TAKU','MACKENZIE','NASS','THJORSA','JOEKULSA A F.',...
'KUSKOKWIM','RHONE','SKEENA','OB','OELFUSA','MEKONG','DANUBE','NELSON RIVER','PO','KAMCHATKA','RHINE','GLOMA','HUANG HE','INDIGIRKA',...
'LULE','RAPEL','SANTA','SKAGIT','KUBAN','TITICACA','NUSHAGAK','BIOBIO','IRRAWADDY','NEGRO','MAJES','CLUTHA','DAULE/VINCES',...
'KALIXAELVEN','MAGDALENA','DRAMSELV','COLVILLE'};
BasinArea=[1139075,1051731,518011,1233148,64959,1024462,829632,28422,49470,423657,51147,30599,...
239678,30760,1745094,258475,668561,191032,5880854,390631,17967,1752001,21211,7527,7311,...
118114,97485,42944,2701040,5678,787256,793704,1099380,73066,54103,190522,42862,988062,341227,...
25127,15689,11882,7961,58935,107215,29513,24108,411516,130062,18612,17118,41993,...
17157,261204,17364,57544];
daysmon=[31,28,31,30,31,30,31,31,30,31,30,31];
load('ModelSetFin.mat') % These are the models and simulations corresponding to the data in ModData,
% these are all models that had the necessary variables not just those in
% Huss and Hock (2018)
InVal=InVal([4,9,16,19,22,24,29,32]); % This limits to just the models that have glacial runoff from Huss and Hock (2018)
% Getting the percent glaciated:
PercGlac=dlmread('GlacialArea.txt')./BasinArea;
GlacArea=dlmread('GlacialArea.txt');
ModNam={'CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GISS-E2-R','INMCM4','MIROC-ESM','NorESM1-M'};


%
%
%
% Loading in the files with glacier runoff
RunoffRcp8p5=nan(length(ModNam),56,1452); 
for kk=1:length(ModNam)
    A=importdata(['discharge_',ModNam{kk},'_rcp85.dat'],' ',2);
    In=A.data;
    for i=1:56
        In2=In(find(In(:,1)==BasinInd(i)),3:14)*10^6;
        RunoffRcp8p5(kk,i,:)=reshape(In2',size(In2,1)*12,1);
    end
end

%
%
%
% Grabbing the data for runoff

% Loading the quality control of the data. Not all of the runoff data is
% good thus did quanitative quality control checks CheckIf has the data
% that passed four check criteria (no nan, reasonable mean value, not zero
% variance, reasonable variance). There are additional basins that did pass
% the quantitative quality control but did not pass the eye test. I set
% these below ("Bad after eye check"). In total 73% of basin model combinations...
load('CheckRunoff.mat')
CheckIf=squeeze(CheckIf(:,:,5));
CheckIf(47,7)=nan; % Bad after eye check
CheckIf(41,1)=nan; % Bad after eye check
CheckIf(30,1)=nan; % Bad after eye check
CheckIf(4,[5 8])=nan; % Bad after eye check
time=[1860:1/12:2101-1/12];
% Rescaling runoff to mm/day and limiting to the time interval
% to 1900 to 2100:
RRcp8p5=DataRcp8p5_RunOff(:,:,find(time==1900):end)*60*60*24; 
time=time(find(time==1900):end);
% Converting from mm/day to m cubed/month: 
% (1/1000) is mm to m;
% repmat(daysmon',[1452/12 1])) is /day to /month;
% (repmat(GlacArea(i),[1452 1])*1000000) is the basin area to get volume water
RRcp8p5Scal=nan(size(RRcp8p5));
for kk=1:length(ModNam)
    for i=1:56
        RRcp8p5Scal(kk,i,:)=squeeze(RRcp8p5(kk,i,:)*(1/1000)).*repmat(daysmon',[2412/12 1]).*(repmat(BasinArea(i),[2412 1])*1000000);
    end
end
RRcp8p5=RRcp8p5Scal;
RRcp8p5wRf=RRcp8p5Scal;

%
%
%
% Adding in the glacial runoff
for i=1:length(ModNam)
    RRcp8p5wRf(i,:,find(time==1980):end)=repmat((1-PercGlac)',[1 1452]).*(squeeze(RRcp8p5wRf(i,:,find(time==1980):end)))+squeeze(RunoffRcp8p5(i,:,:)); % w/ runoff
end

%
%
%
% Calculating the 3 month non parametric RPI for comparison to SPEI:
scalval=[3];
timey=reshape(repmat([1900:2100],[12 1]),12*length([1900:2100]),1);
DataOut=nan(2,length(ModNam),length(BasinNam),length(time));
for i=1:length(ModNam)
    for j=1:length(BasinNam)
            
        % Fixing the output time
        erase_yr=ceil(scalval/12);
        nseas=12;    
            
        % Original Rcp8p5    
        [Z,Z_std]=SPEI_Calc(squeeze(RRcp8p5(i,j,:)),scalval,12,timey,1900,1979);
        DataOut(1,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z; 
        
        % Runoff Rcp8p5
        [Z,Z_std]=SPEI_Calc(squeeze(RRcp8p5wRf(i,j,:)),scalval,12,timey,1900,1979);
        DataOut(2,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z; 
           
    end
end

%

% Doing the comparison:
load('SPEIOut.mat')
tind=[1 612;1081 1692;1801 2412];
CorrVals=nan(2,3,length(ModNam),length(BasinNam));
for i=1:length(ModNam)
    for j=1:length(BasinNam)
        
        if CheckIf(j,i)==1
            
            for h=1:3 % these are the three time periods:

            InOrig=squeeze(SPEIRcp8p5(1,i,j,tind(h,1):tind(h,2)));
            InNew=squeeze(DataOut(1,i,j,tind(h,1):tind(h,2)));
            indOrig=find(isnan(InOrig)==0);
            indNew=find(isnan(InNew)==0);
            ind=intersect(indOrig,indNew);
            if isempty(ind)==0
            CorrVals(1,h,i,j)=corr(InOrig(ind),InNew(ind));
            end
            InOrig=squeeze(SPEIRcp8p5(2,i,j,tind(h,1):tind(h,2)));
            InNew=squeeze(DataOut(2,i,j,tind(h,1):tind(h,2)));
            indOrig=find(isnan(InOrig)==0);
            indNew=find(isnan(InNew)==0);
            ind=intersect(indOrig,indNew);
            if isempty(ind)==0
            CorrVals(2,h,i,j)=corr(InOrig(ind),InNew(ind));
            end
            end
        end

    end
end

%{
% Correlation of future difference:
CorrValsDiff=nan(length(ModNam),length(BasinNam));
for i=1:length(ModNam)
    for j=1:length(BasinNam)
        
        if CheckIf(j,i)==1
            
            for h=3

            InOrig=squeeze(SPEIRcp8p5(1,i,j,tind(h,1):tind(h,2)));
            InNew=squeeze(DataOut(1,i,j,tind(h,1):tind(h,2)));
            InOrig2=squeeze(SPEIRcp8p5(2,i,j,tind(h,1):tind(h,2)));
            InNew2=squeeze(DataOut(2,i,j,tind(h,1):tind(h,2)));
            InOrig=InOrig2-InOrig;
            InNew=InNew2-InNew;
            indOrig=find(isnan(InOrig)==0);
            indNew=find(isnan(InNew)==0);
            ind=intersect(indOrig,indNew);
            if isempty(ind)==0
            CorrValsDiff(i,j)=corr(InOrig(ind),InNew(ind));
            end
            end
        end

    end
end
%}

% Plot the histogram of correlations 
figure(1)
subplot(3,2,[1 2])
set(gca,'LineWidth',2)
InVal=squeeze(CorrVals(1,1,:,:));
histogram(InVal(:),10,'Normalization','probability','FaceColor','k','EdgeColor','w')
hold on
plot([nanmean(InVal(:)) nanmean(InVal(:))],[0 .25],'k--','LineWidth',2)
xlim([-1 1])
ylim([0 .25])
set(gca,'ytick',[.05:.05:.25])
set(gca,'xtick',[-1:.2:1])
set(gca,'yticklabel',{'0.05','0.10','0.15','0.20','0.25'},'FontName','Helvetica','FontSize',18)
title('Histogram of Correlation Coefficients for SRI and SPEI (N=328)','FontName','Helvetica','FontSize',21)
ylabel('Fraction in bin (1900-1950)','FontName','Helvetica','FontSize',21)
subplot(3,2,3)
set(gca,'LineWidth',2)
InVal=squeeze(CorrVals(1,2,:,:));
histogram(InVal(:),10,'Normalization','probability','FaceColor','k','EdgeColor','w')
hold on
plot([nanmean(InVal(:)) nanmean(InVal(:))],[0 .25],'k--','LineWidth',2)
xlim([-1 1])
ylim([0 .25])
set(gca,'ytick',[.05:.05:.25])
set(gca,'xtick',[-1:.2:1])
set(gca,'yticklabel',{'0.05','0.10','0.15','0.20','0.25'},'FontName','Helvetica','FontSize',18)
ylabel('Fraction in bin (1990-2040)','FontName','Helvetica','FontSize',21)
text(-.95,.22,'without glacial runoff','FontName','Helvetica','FontSize',19)
subplot(3,2,4)
set(gca,'LineWidth',2)
InVal=squeeze(CorrVals(2,2,:,:));
histogram(InVal(:),10,'Normalization','probability','FaceColor','k','EdgeColor','w')
hold on
plot([nanmean(InVal(:)) nanmean(InVal(:))],[0 .25],'k--','LineWidth',2)
xlim([-1 1])
ylim([0 .25])
set(gca,'ytick',[.05:.05:.25])
set(gca,'xtick',[-1:.2:1])
set(gca,'yticklabel',{'0.05','0.10','0.15','0.20','0.25'},'FontName','Helvetica','FontSize',18)
text(-.95,.22,'with glacial runoff','FontName','Helvetica','FontSize',19)
subplot(3,2,5)
set(gca,'LineWidth',2)
InVal=squeeze(CorrVals(1,3,:,:));
histogram(InVal(:),10,'Normalization','probability','FaceColor','k','EdgeColor','w')
hold on
plot([nanmean(InVal(:)) nanmean(InVal(:))],[0 .25],'k--','LineWidth',2)
xlim([-1 1])
ylim([0 .25])
set(gca,'ytick',[.05:.05:.25])
set(gca,'xtick',[-1:.2:1])
set(gca,'yticklabel',{'0.05','0.10','0.15','0.20','0.25'},'FontName','Helvetica','FontSize',18)
ylabel('Fraction in bin (2050-2100)','FontName','Helvetica','FontSize',21)
text(-.95,.22,'without glacial runoff','FontName','Helvetica','FontSize',19)
subplot(3,2,6)
set(gca,'LineWidth',2)
InVal=squeeze(CorrVals(2,3,:,:));
histogram(InVal(:),10,'Normalization','probability','FaceColor','k','EdgeColor','w')
hold on
plot([nanmean(InVal(:)) nanmean(InVal(:))],[0 .25],'k--','LineWidth',2)
xlim([-1 1])
ylim([0 .25])
set(gca,'ytick',[.05:.05:.25])
set(gca,'xtick',[-1:.2:1])
set(gca,'yticklabel',{'0.05','0.10','0.15','0.20','0.25'},'FontName','Helvetica','FontSize',18)
text(-.95,.22,'with glacial runoff','FontName','Helvetica','FontSize',19)

% Fixing and saving:

%Enlarging the plot window:
screen_size=get(0,'ScreenSize');
set(figure(1), 'Position',[0 0 screen_size(3) screen_size(4)]);

%Making them bigger:
hAx = findobj(figure(1), 'type', 'axes');
hAx = flipud(hAx);
nF1=1.1;
nF=1.1;
for h = 2:5
    vCurrPos = get(hAx(h), 'position'); % current position
    set(hAx(h), 'position', (vCurrPos.*[1 1 nF1 nF])-[vCurrPos(3)*(nF1-1)/2 vCurrPos(4)*(nF-1)/2  0 0]);
end

%Moving them all down and the left in and the rightin
mind=[0 -.0145 .0145 -.0145 .0145];
for h = 2:5
    vCurrPos = get(hAx(h), 'position'); % current position
    set(hAx(h), 'position', (vCurrPos.*[1 1 1 1])-[mind(h) .07  0 0]);
end

%Making the top taller
hAx = findobj(figure(1), 'type', 'axes');
hAx = flipud(hAx);
nF1=1.005;
nF=1.595;
for h = 1
    vCurrPos = get(hAx(h), 'position'); % current position
    set(hAx(h), 'position', (vCurrPos.*[1 1 nF1 nF])-[vCurrPos(3)*(nF1-1)/2 vCurrPos(4)*(nF-1)/2+.018  0 0]);
end

%
%
%
% Saving 
set(gcf,'color','w');
export_fig Corr_SRI_SPEI -m2.5 -png
export_fig Corr_SRI_SPEI -painters -pdf

function [Out,Out_std]=SPEI_Calc(Data,scale,nseas,yr_vect,y1,y2)

            % Data setting to scaled dataset
            erase_yr=ceil(scale/12);
            A1=[];
            for is=1:scale, A1=[A1,Data(is:length(Data)-scale+is)];end
            XS=sum(A1,2);

            if(scale>1), XS(1:nseas*erase_yr-scale+1)=[];   end
            
            % Trim the year vector
            yr_vect(1:nseas*erase_yr)=[];
            
            % Compute the standardized index
            nn=length(XS);
            SI1=zeros(nn,1);
            for k=1:12
            d=XS(k:12:nn);
            yr_in=yr_vect(k:12:nn);
            syear=find(yr_in==y1);
            if isempty(syear)==1
               syear=1;
            end
            eyear=find(yr_in==y2);
            nnn=length(d);
            bp=zeros(nnn,1);
            for i=1:nnn
            bp(i,1)=sum(d(syear:eyear,1)<=d(i,1));
            end
            y=(bp-0.44)./(length(d(syear:eyear,1))+0.12);
            SI1(k:12:nn,1)=y;
            end
            Out1=norminv(SI1(:,1));
            i_yr=find(yr_vect>=y1 & yr_vect<=y2);
            Out=Out1;
            
            % Setting values below standardization range to -3:
            Out(isnan(Out)==1)=-3;
            
            % Check the standard deviation of the standardization
            % interval:
            Out_std=std(Out(i_yr));
end

