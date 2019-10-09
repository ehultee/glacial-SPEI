%
%
%
% Setin up the basin info
load('InputData.mat') % Total precipitation and variables needed to calculate PET
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
% Getting the percent glaciated:
PercGlac=dlmread('GlacialArea.txt')./BasinArea;
GlacArea=dlmread('GlacialArea.txt');

%
%
%
% Setting up P and PET for the 1980-2100 period

% Calculate PET
PETRcp4p5=nan(length(InVal),length(BasinNam),length(time));
for i=1:length(InVal)
    for j=1:length(BasinNam)
        
        % NETRAD (rsds-rsus+rlds-rlus)
        NETRAD=(squeeze(DataRcp4p5(i,5,j,:))-squeeze(DataRcp4p5(i,6,j,:)))+(squeeze(DataRcp4p5(i,7,j,:))-squeeze(DataRcp4p5(i,8,j,:)));
        
        % VP e=huss*ps/(.622+.378*huss)
        VP=(squeeze(DataRcp4p5(i,4,j,:)).*squeeze(DataRcp4p5(i,2,j,:)))./(0.622*ones(length(time),1)+0.378*squeeze(DataRcp4p5(i,4,j,:)));
        
        % PET
        [PETIn,~]=penmont_model_q(squeeze(DataRcp4p5(i,1,j,:))-273.15,NETRAD,VP,squeeze(DataRcp4p5(i,2,j,:)));
        PETRcp4p5(i,j,:)=PETIn;
        
    end
end
PETRcp8p5=nan(length(InVal),length(BasinNam),length(time));
for i=1:length(InVal)
    for j=1:length(BasinNam)
        
        % NETRAD (rsds-rsus+rlds-rlus)
        NETRAD=(squeeze(DataRcp8p5(i,5,j,:))-squeeze(DataRcp8p5(i,6,j,:)))+(squeeze(DataRcp8p5(i,7,j,:))-squeeze(DataRcp8p5(i,8,j,:)));
        
        % VP e=huss*ps/(.622+.378*huss)
        VP=(squeeze(DataRcp8p5(i,4,j,:)).*squeeze(DataRcp8p5(i,2,j,:)))./(0.622*ones(length(time),1)+0.378*squeeze(DataRcp8p5(i,4,j,:)));
        
        % PET
        [PETIn,~]=penmont_model_q(squeeze(DataRcp8p5(i,1,j,:))-273.15,NETRAD,VP,squeeze(DataRcp8p5(i,2,j,:)));
        PETRcp8p5(i,j,:)=PETIn;
        
    end
end

% Getting P
PRECRcp4p5=squeeze(DataRcp4p5(:,3,:,:))*60*60*24; %kg/m2s to mm/day;
PRECRcp8p5=squeeze(DataRcp8p5(:,3,:,:))*60*60*24; %kg/m2s to mm/day;

% Clipping time on everything to match the glacier outputs
PETRcp4p5=PETRcp4p5(:,:,find(time==1900):end);
PETRcp8p5=PETRcp8p5(:,:,find(time==1900):end);
PRECRcp4p5=PRECRcp4p5(:,:,find(time==1900):end);
PRECRcp8p5=PRECRcp8p5(:,:,find(time==1900):end);
time=single(time(find(time==1900):end));

% Grabbing only the models of interest:
% CanESM2 (4)
% CCSM4 (9)
% CNRM-CM5 (16)
% CSIRO-Mk3-6-0 (19)
% GISS-E2-R (22)
% INMCM4 (24)
% MIROC-ESM (29)
% NorESM1-M (32)
PETRcp4p5=PETRcp4p5([4,9,16,19,22,24,29,32],:,:);
PETRcp8p5=PETRcp8p5([4,9,16,19,22,24,29,32],:,:);
PRECRcp4p5=PRECRcp4p5([4,9,16,19,22,24,29,32],:,:);
PRECRcp8p5=PRECRcp8p5([4,9,16,19,22,24,29,32],:,:);

%
%
%
% Loading in the files with glacier runoff
ModNam={'CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GISS-E2-R','INMCM4','MIROC-ESM','NorESM1-M'};
RunoffRcp4p5=nan(length(ModNam),56,1452);   
for kk=1:length(ModNam)
    A=importdata(['discharge',ModNam{kk},'_rcp45.dat'],' ',2);
    In=A.data;
    for i=1:56
        In2=In(find(In(:,1)==BasinInd(i)),3:14)*10^6;
        RunoffRcp4p5(kk,i,:)=reshape(In2',size(In2,1)*12,1);
    end
end
RunoffRcp8p5=nan(length(ModNam),56,1452); 
for kk=1:length(ModNam)
    A=importdata(['discharge',ModNam{kk},'_rcp85.dat'],' ',2);
    In=A.data;
    for i=1:56
        In2=In(find(In(:,1)==BasinInd(i)),3:14)*10^6;
        RunoffRcp8p5(kk,i,:)=reshape(In2',size(In2,1)*12,1);
    end
end

% Scaling to correct units <-- Make sure the area part is correct, % dividing by the basin area to go from volume, 1000000 is to get to meters2 from km2
for kk=1:length(ModNam)
    for i=1:56
    RunoffRcp4p5(kk,i,:)=squeeze(RunoffRcp4p5(kk,i,:))./(repmat(GlacArea(i),[1452 1])*1000000); % dividing by the basin area to go from volume, 1000000 is to get to meters2 from km2
    RunoffRcp8p5(kk,i,:)=squeeze(RunoffRcp8p5(kk,i,:))./(repmat(GlacArea(i),[1452 1])*1000000);  % dividing by the basin area to go from volume, 1000000 is to get to meters2 from km2
    end
end

% getting things to mm/day from m/month
for kk=1:length(ModNam)
    for i=1:56
    RunoffRcp4p5(kk,i,:)=(squeeze(RunoffRcp4p5(kk,i,:))./repmat(daysmon',[1452/12 1]))*1000; 
    RunoffRcp8p5(kk,i,:)=(squeeze(RunoffRcp8p5(kk,i,:))./repmat(daysmon',[1452/12 1]))*1000; 
    end
end

%
%
%
% Setting up the moisture balance
DRcp4p5=PRECRcp4p5-PETRcp4p5;
DRcp8p5=PRECRcp8p5-PETRcp8p5;
DRcp4p5wRf=PRECRcp4p5-PETRcp4p5;
DRcp8p5wRf=PRECRcp8p5-PETRcp8p5;

% Have to fix this...
for i=1:length(ModNam)
DRcp4p5wRf(i,:,find(time==1980):end)=repmat((1-PercGlac)',[1 1452]).*(squeeze(PRECRcp4p5(i,:,find(time==1980):end))-squeeze(PETRcp4p5(i,:,find(time==1980):end)))+repmat(PercGlac',[1 1452]).*squeeze(RunoffRcp4p5(i,:,:)); % w/ runoff
DRcp8p5wRf(i,:,find(time==1980):end)=repmat((1-PercGlac)',[1 1452]).*(squeeze(PRECRcp8p5(i,:,find(time==1980):end))-squeeze(PETRcp8p5(i,:,find(time==1980):end)))+repmat(PercGlac',[1 1452]).*squeeze(RunoffRcp8p5(i,:,:)); % w/ runoff
end
% So this is percentage glacier time runoff but percentage not glacier time
% precipitation minus the PET


%
%
%
% Calculating SPEI 3,7,11,15,19,23,27 months:
scalvals=[3,7,11,15,19,23,27];
timey=reshape(repmat([1900:2100],[12 1]),12*length([1900:2100]),1);
SPEIRcp4p5=nan(2,length(scalvals),length(ModNam),length(BasinNam),length(time));
SPEIRcp8p5=nan(2,length(scalvals),length(ModNam),length(BasinNam),length(time));
for i=1:length(ModNam)
    for j=1:length(BasinNam)
        for k=1:length(scalvals)
        
        % Fixing the output time
        erase_yr=ceil(scalvals(k)/12);
        nseas=12;
        
        % Original Rcp4p5    
        [Z,Z_std]=SPEI_time3(squeeze(DRcp4p5(i,j,:)),scalvals(k),12,timey,1900,1979);
        SPEIRcp4p5(1,k,i,j,nseas*erase_yr-scalvals(k)+1:end-scalvals(k))=Z;
        Z_std % Should always be ~1
        
        % Runoff Rcp4p5
        [Z,Z_std]=SPEI_time3(squeeze(DRcp4p5wRf(i,j,:)),scalvals(k),12,timey,1900,1979);
        SPEIRcp4p5(2,k,i,j,nseas*erase_yr-scalvals(k)+1:end-scalvals(k))=Z;
        Z_std
        
        % Original Rcp8p5
        [Z,Z_std]=SPEI_time3(squeeze(DRcp8p5(i,j,:)),scalvals(k),12,timey,1900,1979);
        SPEIRcp8p5(1,k,i,j,nseas*erase_yr-scalvals(k)+1:end-scalvals(k))=Z;
        Z_std
        
        % Runoff Rcp8p5
        [Z,Z_std]=SPEI_time3(squeeze(DRcp8p5wRf(i,j,:)),scalvals(k),12,timey,1900,1979);
        SPEIRcp8p5(2,k,i,j,nseas*erase_yr-scalvals(k)+1:end-scalvals(k))=Z;
        Z_std
        
        end
    end
end

%
%
%
% Saving the output:
save SPEIOut SPEIRcp4p5 SPEIRcp8p5 timey time scalvals ModNam
