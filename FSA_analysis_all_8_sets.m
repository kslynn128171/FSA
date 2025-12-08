clearvars
clc
close all
fclose all;
opt.TINY=0.5; % ratio threshold for "tiny" peaks
opt.mass=0;
opt.tol=0.003;
ppm=[0.0013,0.0006,0.0005,0.0029,0.0013,0.0011,0.0011,0.0005];
opt.MaxPeak=30;
outputname='Supplementary_File1_8_sets.xlsx';
datasets={'HMDB_QTOF_spectra';
    'HMDB_ITFT_spectra';
    'HMDB_QFT_spectra';
    'GNPS_NPL_MAXIS_spectra';
    'GNPS_NPL_QTOF_spectra';
    'GNPS_NPL_ORBITRAP_spectra';
    'PNNLLIPID_spectra';
    'CASMI2016_cat2and3'};
sheetnames={'HMDB_QTOF','HMDB_ITFT','HMDB_QFT','NPL_MAXIS','NPL_QTOF','NPL_OT','LIPID_OT','CASMI_OT'};
noelem=[0,0,0,6,6,6,5,0];
varNames = {'No','Name','Instrument_Type','Ionization_Mode','Mass','Formula','Rank_Reduced','Original_Rank','Original_Candidates','Peak_Num',...
        'Frag_Num','Frag_w_Formula','Score','Is_Reduced','Is_Top1_Uniq','Char_Frag_Percent','Is_CFP_Max','Is_CFP_Max_Uniq','Top2_Sore_Diff',...
        'Peak_Dist_SD','Comp_Time','Note'};
varTypes = {'double','string','string','string','double','string','double','double','double','double',...
        'double','double','double','double','double','double','double','double','double','double',...
        'double','string'};
dlen=length(datasets);
% prepare summary of FSM results
vName={};
for i=1:dlen
    vName=[vName,{['SpectrumNum',num2str(i),'(Percentage)'],['CompoundNum',num2str(i),'(Percentage)']}];
end
vType={'string','string'};
vType=repmat(vType,1,dlen);
rName={'Total','Problematic','No_Answer_within_Candidates','No_MDR','Solo_Candidate','Solo_Ran#k1','Top-3_Rank','Top-5_Rank'};
rank1 = table('Size',[length(rName) length(vName)],'VariableTypes',vType,'VariableNames',vName,'RowNames',rName);
totaltime=0;
spectvalues=zeros(dlen,5);
metabvalues=zeros(dlen,5);
for d=1:dlen
    starttime=cputime;
    load PubChemMetabolite.mat;
    % load test spectra
    load(datasets{d});
    mmw=extractfield(spectra,'monisotopic_molecular_weight');
    if ~isnumeric(mmw)
        maxmass=max(str2double(mmw)); % max molecular weight in the spectra
    else
        maxmass=max(mmw); % max molecular weight in the spectra
    end
    if noelem(d)>0
        useid=(mass <= (maxmass+1.1)) & (elemnum(:,noelem(d))==0); % remove compounds with P
    else
        useid=(mass <= (maxmass+1.1));
    end
    elemnum=uint8(elemnum(useid,:));
    formula=formula(useid);
    mass=mass(useid);
    % Build parallel Constant from data
    spectnum=size(spectra,1); % total number of spectra
    opt.MatchTolerance=ppm(d);% mass match tolerance 
    spect = table('Size',[spectnum length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    usedtime=cputime-starttime;
    disp(['Dataset ',num2str(d),' used ',num2str(usedtime),' secs for database loading.'])
    for i=1:spectnum % for each spectrum
        if d <=3
            msg=SingleSpectrumIdent_HMDB_massmatch(i,elemnum,formula,mass,spectra(i),opt);
        elseif d < 8
            msg=SingleSpectrumIdent_NatProd_massmatch(i,elemnum,formula,mass,spectra(i),opt);
        else
            msg=SingleSpectrumIdent_CASMI2016_massmatch(i,elemnum,formula,mass,spectra(i),opt);
        end
        spect(i,:)=msg;
    end
    usedtime=cputime-starttime;
    totaltime=totaltime+usedtime;
    disp(['Dataset ',datasets{d},' used ',num2str(usedtime),' secs.'])
    writetable(spect,outputname,'Sheet',[sheetnames{d},'_result']);
    Rank_Reduced=spect.Rank_Reduced;
    Original_Rank=spect.Original_Rank;
    Original_Candidates=spect.Original_Candidates;
    Frag_w_Formula=spect.Frag_w_Formula;
    Char_Frag_Percent=spect.Char_Frag_Percent;
    Score=spect.Score;
    % problematic
    prob=(Rank_Reduced==-1) & (Original_Rank==-1) & (Original_Candidates==-1);
    % no formula matches with the correct formula
    nomatch=(Rank_Reduced==-1) & (Original_Rank==-1) & (Original_Candidates>0);
    % unidentifiable
    unident=prob | nomatch;
    % no MDR
    noMDall=(Rank_Reduced>=1) & (Original_Rank>=1) & (Original_Candidates>1) & (Char_Frag_Percent==-1);
    % solo and correct candidate
    solosol=(Rank_Reduced==1) & (Original_Rank==1) & (Original_Candidates==1);
    inapplicable=prob | noMDall | solosol;
    % unique rank1 is the correct formula
    solorank1=(Rank_Reduced==1) & (Original_Candidates>=1);
    top3sol=(Rank_Reduced<=3) & (Rank_Reduced>0);
    top5sol=(Rank_Reduced<=5) & (Rank_Reduced>0);
    totalspect=size(spect,1);
    totalcomp=size(unique(spect(:,2)),1);
    inapplicablecomp=size(unique(spect(inapplicable,2)),1);
    probspect=sum(prob);
    probcomp=size(unique(spect(prob,2)),1);
    nomatchspect=sum(nomatch);
    nomatchcomp=size(unique(spect(nomatch,2)),1);
    noMDspect=sum(noMDall);
    noMDcomp=size(unique(spect(noMDall,2)),1);
    SoloSolspect=sum(solosol);
    SoloSolcomp=size(unique(spect(solosol,2)),1);
    is_applicable=~inapplicable;
    applicablespect=sum(is_applicable);
    applicablecomp=size(unique(spect(is_applicable,2)),1);
    solorank1spect=sum(solorank1);
    solorank1comp=size(unique(spect(solorank1,2)),1);
    top3solspect=sum(top3sol);
    top3solcomp=size(unique(spect(top3sol,2)),1);
    top5solspect=sum(top5sol);
    top5solcomp=size(unique(spect(top5sol,2)),1);
    rank1(1,(2*(d-1)+1):(2*d))={num2str(totalspect),num2str(totalcomp)};
    rank1(2,(2*(d-1)+1):(2*d))={num2str(probspect),num2str(probcomp)};
    rank1(3,(2*(d-1)+1):(2*d))={num2str(nomatchspect),num2str(nomatchcomp)};
    rank1(4,(2*(d-1)+1):(2*d))={num2str(noMDspect),num2str(noMDcomp)};
    rank1(5,(2*(d-1)+1):(2*d))={[num2str(SoloSolspect),'(',num2str(SoloSolspect/totalspect*100,'%6.2f'),'%)'],...
        [num2str(SoloSolcomp),'(',num2str(SoloSolcomp/totalcomp*100,'%6.2f'),'%)']};
    rank1(6,(2*(d-1)+1):(2*d))={[num2str(solorank1spect),'(',num2str(solorank1spect/totalspect*100,'%6.2f'),'%)'],...
        [num2str(solorank1comp),'(',num2str(solorank1comp/totalcomp*100,'%6.2f'),'%)']};
    rank1(7,(2*(d-1)+1):(2*d))={[num2str(top3solspect),'(',num2str(top3solspect/totalspect*100,'%6.2f'),'%)'],...
        [num2str(top3solcomp),'(',num2str(top3solcomp/totalcomp*100,'%6.2f'),'%)']};
    rank1(8,(2*(d-1)+1):(2*d))={[num2str(top5solspect),'(',num2str(top5solspect/totalspect*100,'%6.2f'),'%)'],...
        [num2str(top5solcomp),'(',num2str(top5solcomp/totalcomp*100,'%6.2f'),'%)']};
    spectvalues(d,1)=totalspect;
    spectvalues(d,2)=SoloSolspect;
    spectvalues(d,3)=solorank1spect;
    spectvalues(d,4)=top3solspect;
    spectvalues(d,5)=top5solspect;
    metabvalues(d,1)=totalcomp;
    metabvalues(d,2)=SoloSolcomp;
    metabvalues(d,3)=solorank1comp;
    metabvalues(d,4)=top3solcomp;
    metabvalues(d,5)=top5solcomp;
end
writetable(rank1,outputname,'Sheet','ident_summary','WriteRowNames',true);
totaltime
ss=sum(spectvalues,1);
mm=sum(metabvalues,1);
disp(ss);
disp(num2str(ss(2:end)/ss(1)*100,',%5.2f'))
disp(mm);
disp(num2str(mm(2:end)/mm(1)*100,',%5.2f'))