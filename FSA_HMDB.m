clearvars
clc
close all
fclose all;
% ---------------------------------------------------------------------
% This program is developed for the generation of Table 3 and the
% Supplementary File 1 for the paper entiled "Automated Metabolite 
% Candidate Ranking Using Formula Subset Analysis for Characteristic 
% Fragments in Liquid Chromatography–Tandem Mass Spectrometry".
%
% Author: Ke-Shiuan Lynn Ph.D.
% Assistant Professor
% Department of Mathematics
% Fu-Jen Catholic University
% Email: 128171@mail.fju.edu.tw
% Final Update: Nov. 12, 2021
% ---------------------------------------------------------------------
% Parameter settings
opt.TINY=0; % ratio threshold for "tiny" peaks
opt.PPM=0.001;% mass match tolerance
% create a progress bar
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
% load test spectra
load PubChemMetabolite.mat;
output_filename='FSA_HMDB_result_CHONSP.xlsx'; % file name for the identification result
HMDBTable=readtable('HMDB_CHONPS_MSMS.xlsx','FileType','spreadsheet',...
   'TextType','string','VariableNamingRule','modify');
% remove Pubmed data whose molecular weight is larger than the max molecular weight in the spectra
if ~isnumeric(HMDBTable.monisotopic_molecular_weight)
    maxmass=max(str2double(HMDBTable.monisotopic_molecular_weight)); % max molecular weight in the spectra
else
    maxmass=max(HMDBTable.monisotopic_molecular_weight); % max molecular weight in the spectra
end
useid=mass <= (maxmass+1.1);
elemnum=uint8(elemnum(useid,:));
formula=formula(useid);
mass=mass(useid);
spectra=table2struct(HMDBTable);
spectnum=size(HMDBTable,1); % total number of spectra
% construct a table to record the identification information
spect = table('Size',[spectnum length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
idlen=size(spect,1);
% initialize the progress bar
h = waitbar(0, 'Waiting...','tag','waitfig');
% start the identification for each MS/MS spectrum
starttime=cputime;
for i=1:idlen % for each spectrum
    msg=SingleSpectrumIdent_submit(i,elemnum,formula,mass,spectra(i),opt);
    spect(i,:)=msg;
    ratio=1.0*i/idlen;
    waitbar(ratio, h, ['In HMDB: ',num2str(ratio*100,'%6.2f'),'% (',num2str(i),'/',num2str(idlen),') finished']);
end
totaltime=cputime-starttime
% write FSM results
writetable(spect,output_filename,'Sheet','HMDB_result');
% prepare summary of FSM results
vName={'SpectrumNum(Percentage)','CompoundNum(Percentage)'};
vType={'string','string'};
rName={'Total','Problematic','No_Answer_within_Candidates','No_MDR','Total_Identifiable','Solo_Candidate','Solo_Ran#k1','Top-3_Rank','Top-5_Rank'};
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
writetable(rank1,output_filename,'Sheet','HMDB_summary','WriteRowNames',true);
%save FSA_HMDB_result_all_tables_001_new.mat spect rank1;
save FSA_HMDB_result_CHONSP_Halogen_001.mat spect rank1;
delete(h);
