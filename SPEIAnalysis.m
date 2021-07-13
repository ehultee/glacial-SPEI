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
InVal=InVal([4,9,16,19,22,24,29,32]); % This limits to just the models that have glacial runoff from Huss and Hock (2018)
% Getting the percent glaciated:
PercGlac=dlmread('GlacialArea.txt')./BasinArea;
GlacArea=dlmread('GlacialArea.txt');
ModNam={'CanESM2','CCSM4','CNRM-CM5','CSIRO-Mk3-6-0','GISS-E2-R','INMCM4','MIROC-ESM','NorESM1-M'};

%
%
%
% Reading in and processing the global carbon dioxide emissions
Rcp45CDIO=xlsread('RCP45_CDIO.xls');
Rcp85CDIO=xlsread('RCP85_CDIO.xls');
Rcp45year=Rcp45CDIO(:,1)';
Rcp85year=Rcp85CDIO(:,1)';
Rcp45CDIO=Rcp45CDIO(find(Rcp45year==1860):find(Rcp45year==2100),2);
Rcp85CDIO=Rcp85CDIO(find(Rcp45year==1860):find(Rcp45year==2100),2);
Rcp45CDIO=repmat(Rcp45CDIO,1,12);
Rcp45CDIO=reshape(Rcp45CDIO',2892,1);
Rcp85CDIO=repmat(Rcp85CDIO,1,12);
Rcp85CDIO=reshape(Rcp85CDIO',2892,1);

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
        [PETIn,~]=penmont_model_q_conduct(squeeze(DataRcp4p5(i,1,j,:))-273.15,NETRAD,VP,squeeze(DataRcp4p5(i,2,j,:)),Rcp45CDIO);
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
        [PETIn,~]=penmont_model_q_conduct(squeeze(DataRcp8p5(i,1,j,:))-273.15,NETRAD,VP,squeeze(DataRcp8p5(i,2,j,:)),Rcp85CDIO);
        PETRcp8p5(i,j,:)=PETIn;
        
    end
end

% Rescaling P to the correct units and saving to a specific variable
PRECRcp4p5=squeeze(DataRcp4p5(:,3,:,:))*60*60*24; %kg/m2s to mm/day;
PRECRcp8p5=squeeze(DataRcp8p5(:,3,:,:))*60*60*24; %kg/m2s to mm/day;

% Clipping time on everything to match the glacier outputs
PETRcp4p5=PETRcp4p5(:,:,find(time==1900):end);
PETRcp8p5=PETRcp8p5(:,:,find(time==1900):end);
PRECRcp4p5=PRECRcp4p5(:,:,find(time==1900):end);
PRECRcp8p5=PRECRcp8p5(:,:,find(time==1900):end);
time=single(time(find(time==1900):end));

%
%
%
% Converting from mm/day to m cubed/month: 
% (1/1000) is mm to m;
% repmat(daysmon',[1452/12 1])) is /day to /month;
% (repmat(GlacArea(i),[1452 1])*1000000) is the basin area to get volume water
PETRcp4p5Scal=nan(size(PETRcp4p5));
PETRcp8p5Scal=nan(size(PETRcp8p5));
PRECRcp4p5Scal=nan(size(PRECRcp4p5));
PRECRcp8p5Scal=nan(size(PRECRcp8p5));
for kk=1:length(ModNam)
    for i=1:56
        PETRcp4p5Scal(kk,i,:)=squeeze(PETRcp4p5(kk,i,:)*(1/1000)).*repmat(daysmon',[2412/12 1]).*(repmat(BasinArea(i),[2412 1])*1000000);
        PETRcp8p5Scal(kk,i,:)=squeeze(PETRcp8p5(kk,i,:)*(1/1000)).*repmat(daysmon',[2412/12 1]).*(repmat(BasinArea(i),[2412 1])*1000000);
        PRECRcp4p5Scal(kk,i,:)=squeeze(PRECRcp4p5(kk,i,:)*(1/1000)).*repmat(daysmon',[2412/12 1]).*(repmat(BasinArea(i),[2412 1])*1000000);
        PRECRcp8p5Scal(kk,i,:)=squeeze(PRECRcp8p5(kk,i,:)*(1/1000)).*repmat(daysmon',[2412/12 1]).*(repmat(BasinArea(i),[2412 1])*1000000);
    end
end
PETRcp4p5=PETRcp4p5Scal;
PETRcp8p5=PETRcp8p5Scal;
PRECRcp4p5=PRECRcp4p5Scal;
PRECRcp8p5=PRECRcp8p5Scal;

%
%
%
% Loading in the files with glacier runoff
RunoffRcp4p5=nan(length(ModNam),56,1452);   
for kk=1:length(ModNam)
    A=importdata(['discharge_',ModNam{kk},'_rcp45.dat'],' ',2);
    In=A.data;
    for i=1:56
        In2=In(find(In(:,1)==BasinInd(i)),3:14)*10^6;
        RunoffRcp4p5(kk,i,:)=reshape(In2',size(In2,1)*12,1);
    end
end
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
% Setting up the moisture balance
DRcp4p5=PRECRcp4p5-PETRcp4p5;
DRcp8p5=PRECRcp8p5-PETRcp8p5;
DRcp4p5wRf=PRECRcp4p5-PETRcp4p5;
DRcp8p5wRf=PRECRcp8p5-PETRcp8p5;

%
%
%
% Adding in the glacial runoff
for i=1:length(ModNam)
    DRcp4p5wRf(i,:,find(time==1980):end)=repmat((1-PercGlac)',[1 1452]).*(squeeze(PRECRcp4p5(i,:,find(time==1980):end)))+squeeze(RunoffRcp4p5(i,:,:))-squeeze(PETRcp4p5(i,:,find(time==1980):end)); % w/ runoff
    DRcp8p5wRf(i,:,find(time==1980):end)=repmat((1-PercGlac)',[1 1452]).*(squeeze(PRECRcp8p5(i,:,find(time==1980):end)))+squeeze(RunoffRcp8p5(i,:,:))-squeeze(PETRcp8p5(i,:,find(time==1980):end)); % w/ runoff
end

%
%
%
% Calculating the 3 month non parametric SPEI:
scalval=3;
timey=reshape(repmat([1900:2100],[12 1]),12*length([1900:2100]),1);
SPEIRcp4p5=nan(2,length(ModNam),length(BasinNam),length(time));
SPEIRcp8p5=nan(2,length(ModNam),length(BasinNam),length(time));
for i=1:length(ModNam)
    for j=1:length(BasinNam)
        
        % Fixing the output time
        erase_yr=ceil(scalval/12);
        nseas=12;
        
        % Original Rcp4p5    
        [Z,Z_std]=SPEI_Calc(squeeze(DRcp4p5(i,j,:)),scalval,12,timey,1900,1979);
        SPEIRcp4p5(1,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z;
        Z_std % Should always be ~1
        
        % Runoff Rcp4p5
        [Z,Z_std]=SPEI_Calc(squeeze(DRcp4p5wRf(i,j,:)),scalval,12,timey,1900,1979);
        SPEIRcp4p5(2,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z;
        Z_std
        
        % Original Rcp8p5
        [Z,Z_std]=SPEI_Calc(squeeze(DRcp8p5(i,j,:)),scalval,12,timey,1900,1979);
        SPEIRcp8p5(1,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z;
        Z_std
        
        % Runoff Rcp8p5
        [Z,Z_std]=SPEI_Calc(squeeze(DRcp8p5wRf(i,j,:)),scalval,12,timey,1900,1979);
        SPEIRcp8p5(2,i,j,nseas*erase_yr-scalval+1:end-scalval)=Z;
        Z_std
        
    end
end

%
%
%
% Saving the output:
save SPEIOut SPEIRcp4p5 SPEIRcp8p5 timey time  ModNam

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
