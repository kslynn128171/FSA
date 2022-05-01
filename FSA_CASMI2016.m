clearvars
clc
close all
opt.TINY=0; % ratio threshold for "tiny" peaks
opt.PPM=0.001;% mass match tolerance for QTOF, ITFT, and QFT
% load PubChem compound database for mass and formula
load PubChemMetabolite.mat;
solution_name='solutions_casmi2016_cat2and3_CHONSP.csv';
output_filename='FSA_result_CASMI2016_Cat2and3_Challenge_CHONSP.xlsx'; % file name for the identification result
%load PubChemMetabolite_add_halogen.mat
%solution_name='CASMI2016_Cat2and3_Train_Valid_all.xlsx';
%output_filename='FSA_result_CASMI2016_Cat2and3_Challenge_CHONSP_Halogen.xlsx'; % file name for the identification result
peakdata_dirname='CASMI2016_Cat2and3_All';

hdl = findobj('tag','waitfig');
if ~isempty(hdl) %delete the old progress bar if exist
    delete(hdl);
end
% The output file columns and their variable types:
varNames = {'No','Name','Instrument_Type','Ionization_Mode','Mass','Formula','Rank_Reduced','Original_Rank','Original_Candidates','Peak_Num',...
        'Frag_Num','Frag_w_Formula','Score','Is_Reduced','Is_Top1_Uniq','Char_Frag_Percent','Is_CFP_Max','Is_CFP_Max_Uniq','Top2_Sore_Diff',...
        'Peak_Dist_SD','Comp_Time','Note'};
varTypes = {'double','string','string','string','double','string','double','double','double','double',...
        'double','double','double','double','double','double','double','double','double','double',...
        'double','string'};
Tb = readtable(solution_name,...
    'FileType','spreadsheet','TextType','string','VariableNamingRule','modify',...
    'VariableNamingRule','preserve');
% remove Pubmed data whose molecular weight is larger than the max molecular weight in the spectra
maxmass=max(Tb.PRECURSOR_MZ); % max molecular weight in the spectra
useid=mass <= (maxmass+1.1);
elemnum=uint8(elemnum(useid,:));
formula=formula(useid);
mass=mass(useid);
% remove records with percursor types other than [M+H]+ and [M-H]-
isProtonated=strcmp(Tb.ION_MODE,'POSITIVE');
spectnum=size(Tb,1);
Tb.ms2_peak=cell(spectnum,1);
Tb.instrument_type=cell(spectnum,1);
Tb.monisotopic_molecular_weight=cell(spectnum,1);
Tb.SMILES=[];
%Tb.INCHI=[];
%Tb.INCHIKEY=[];
%Tb.IUPAC=[];
%Tb.CSID=[];
%Tb.PC_CID=[];
% rearrange spectrum data 
for i=1:spectnum
    if strcmpi(Tb.ION_MODE{i},'positive')
        Tb.monisotopic_molecular_weight{i}=num2str(Tb.PRECURSOR_MZ(i)-1.007276);
    else
        Tb.monisotopic_molecular_weight{i}=num2str(Tb.PRECURSOR_MZ(i)+1.007276);
    end
    Tb.instrument_type{i}='unknown';
    tempmtx=readmatrix([peakdata_dirname,'\',lower(Tb.ION_MODE{i}),'\',Tb.ChallengeName{i},'.txt'],'Delimiter','\t');
    ms2mtx=[zeros(Tb.nPeaks(i),1),tempmtx(:,[2,1]),zeros(Tb.nPeaks(i),1)];
    ms2vec=reshape(ms2mtx',1,numel(ms2mtx));
    ms2str=regexprep(num2str(ms2vec),'\s+',';');
    Tb.ms2_peak{i}=ms2str;
end
Tb = renamevars(Tb,["FORMULA","NAME","ION_MODE"],["chemical_formula","name","ionization_mode"]);
spectra=table2struct(Tb);
% construct a table to record the identification information
spect = table('Size',[spectnum length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
% initialize the progress bar
h = waitbar(0, 'Waiting...','tag','waitfig');
% start the identification for each MS/MS spectrum
starttime=cputime;
for i=1:spectnum % for each spectrum
    msg = SingleSpectrumIdent_submit(i,elemnum,formula,mass,spectra(i),opt);
    spect(i,:)=msg;
    ratio=1.0*i/spectnum;
    waitbar(ratio, h, ['In CASMI: ',num2str(ratio*100,'%6.2f'),'% (',num2str(i),'/',num2str(spectnum),') finished']);
end
totaltime=cputime-starttime
% write FSA results
writetable(spect,output_filename,'Sheet','ID_result');
% prepare summary of FSA results
vName={'SpectrumNum(Percentage)','CompoundNum(Percentage)'};
vType={'string','string'};
rName={'Total','Problematic','No_Answer_within_Candidates','No_MDR','Total_Identifiable','Solo_Candidate','Solo_Ran#k1','Total_Rank1','Top-5_Rank'};
rank1 = table('Size',[length(rName) length(vName)],'VariableTypes',vType,'VariableNames',vName,'RowNames',rName);
Rank_Reduced=spect.Rank_Reduced;
Original_Rank=spect.Original_Rank;
Original_Candidates=spect.Original_Candidates;
Frag_w_Formula=spect.Frag_w_Formula;
Score=spect.Score;
prob=(Rank_Reduced==-1) & (Original_Rank==-1) & (Original_Candidates==-1);
nomatch=(Rank_Reduced==-1) & (Original_Rank==-1) & (Original_Candidates>0) & (Score==-1) & (Frag_w_Formula == -1);
noMDall=(Rank_Reduced==-1) & (Original_Rank==-1) & (Original_Candidates>0) & (Score==-1) & (Frag_w_Formula == 0);
unident=prob | nomatch | noMDall;
solosol=(Rank_Reduced==1) & (Original_Rank==1) & (Original_Candidates==1);
solorank1=(Rank_Reduced==1) & (Original_Candidates>=1);
%allrank1=(Rank_Reduced==1);
top3sol=(Rank_Reduced<=3) & (Rank_Reduced>0);
top5sol=(Rank_Reduced<=5) & (Rank_Reduced>0);
totalspect=size(spect,1);
totalcomp=size(unique(spect(:,2)),1);
unidentcomp=size(unique(spect(unident,2)),1);
probspect=sum(prob);
probcomp=size(unique(spect(prob,2)),1);
nomatchspect=sum(nomatch);
nomatchcomp=size(unique(spect(nomatch,2)),1);
noMDspect=sum(noMDall);
noMDcomp=size(unique(spect(noMDall,2)),1);
SoloSolspect=sum(solosol);
SoloSolcomp=size(unique(spect(solosol,2)),1);
is_ident=~unident;
identspect=sum(is_ident);
identcomp=size(unique(spect(is_ident,2)),1);
solorank1spect=sum(solorank1);
solorank1comp=size(unique(spect(solorank1,2)),1);
top3solspect=sum(top3sol);
top3solcomp=size(unique(spect(top3sol,2)),1);
top5solspect=sum(top5sol);
top5solcomp=size(unique(spect(top5sol,2)),1);
rank1(1,:)={num2str(totalspect),num2str(totalcomp)};
rank1(2,:)={num2str(probspect),num2str(probcomp)};
rank1(3,:)={num2str(nomatchspect),num2str(nomatchcomp)};
rank1(4,:)={num2str(noMDspect),num2str(noMDcomp)};
rank1(5,:)={num2str(identspect),num2str(identcomp)};
rank1(6,:)={[num2str(SoloSolspect),'(',num2str(SoloSolspect/identspect*100,'%6.2f'),'%)'],...
    [num2str(SoloSolcomp),'(',num2str(SoloSolcomp/identcomp*100,'%6.2f'),'%)']};
rank1(7,:)={[num2str(solorank1spect),'(',num2str(solorank1spect/identspect*100,'%6.2f'),'%)'],...
    [num2str(solorank1comp),'(',num2str(solorank1comp/identcomp*100,'%6.2f'),'%)']};
rank1(8,:)={[num2str(top3solspect),'(',num2str(top3solspect/identspect*100,'%6.2f'),'%)'],...
    [num2str(top3solcomp),'(',num2str(top3solcomp/identcomp*100,'%6.2f'),'%)']};
rank1(9,:)={[num2str(top5solspect),'(',num2str(top5solspect/identspect*100,'%6.2f'),'%)'],...
    [num2str(top5solcomp),'(',num2str(top5solcomp/identcomp*100,'%6.2f'),'%)']};
% write summary
writetable(rank1,output_filename,'Sheet','summary','WriteRowNames',true);
save FSA_CASMI2016_CHONSP_Halogen_result_challenge_tables_001.mat spect rank1
%save FSA_CASMI2016_CHONSP_result_All_tables_001.mat spect rank1
%save FSA_CASMI2016_CHONSP_Halogen_result_All_tables_001.mat spect rank1
delete(h);