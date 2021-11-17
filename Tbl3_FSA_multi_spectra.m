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
output_filename='FSA_result_all_spect.xlsx'; % file name for the identification result
opt.TINY=0; % ratio threshold for "tiny" peaks

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
del=[0.0012 0.0011 0.0008];% mass match tolerance for QTOF, ITFT, and QFT
% the test datasets
datasets={'FSA_test_QTOF','FSA_test_ITFT','FSA_test_QFT'};
instrument={'QTOF','ITFT','QFT'};
for d=1:3
    % load test spectra
    eval(['load ',datasets{d}]);
    opt.PPM=del(d);% mass match tolerance
    % remove Pubmed data whose molecular weight is larger than the max molecular weight in the spectra
    maxmass=max(str2double(test1.monisotopic_molecular_weight)); % max molecular weight in the spectra
    useid=mass <= (maxmass+1.1);
    elemnum=uint8(elemnum(useid,:));
    formula=formula(useid);
    mass=mass(useid);
    spectra=table2struct(test1);
    spectnum=size(test1,1); % total number of spectra
    % construct a table to record the identification information
    spect = table('Size',[spectnum length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
    idlen=size(spect,1);
    % initialize the progress bar
    h = waitbar(0, 'Waiting...','tag','waitfig');
    % start the identification for each MS/MS spectrum
    for i=1:idlen % for each spectrum
        msg=SingleSpectrumIdent(i,elemnum,formula,mass,spectra(i),opt);
        spect(i,:)=msg;
        ratio=1.0*i/idlen;
        waitbar(ratio, h, ['In ',instrument{d},': ',num2str(ratio*100,'%6.2f'),'% (',num2str(i),'/',num2str(idlen),') finished']);
    end
    % write FSM results
    writetable(spect,output_filename,'Sheet',[instrument{d},'result']);
    % prepare summary of FSM results
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
    solorank1=(Rank_Reduced==1) & (Original_Rank==1) & (Original_Candidates>=1);
    allrank1=(Rank_Reduced==1);
    top5sol=(Rank_Reduced<=5) & (Rank_Reduced>0) & (Original_Rank<=5);
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
    allrank1spect=sum(allrank1);
    allrank1comp=size(unique(spect(allrank1,2)),1);
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
    rank1(8,:)={[num2str(allrank1spect),'(',num2str(allrank1spect/identspect*100,'%6.2f'),'%)'],...
        [num2str(allrank1comp),'(',num2str(allrank1comp/identcomp*100,'%6.2f'),'%)']};
    rank1(9,:)={[num2str(top5solspect),'(',num2str(top5solspect/identspect*100,'%6.2f'),'%)'],...
        [num2str(top5solcomp),'(',num2str(top5solcomp/identcomp*100,'%6.2f'),'%)']};
    % write summary
    writetable(rank1,output_filename,'Sheet',[instrument{d},'_summary'],'WriteRowNames',true);
    eval(['save FSA_',instrument{d},'_result_tables.mat spect rank1']);
    delete(h);
end
