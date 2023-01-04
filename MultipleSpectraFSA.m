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
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please fill in the following required information
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter settings
opt.TINY=0.5; % ratio threshold for "tiny" peaks
opt.PPM=0.001;% mass match tolerance
opt.MaxRankNumber=1;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please provide your peak list and output filename
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MS2Table=readtable('MSMS_peaklist_example.xlsx','FileType','spreadsheet',...
   'TextType','string','VariableNamingRule','modify'); % file name of the spectral peak list
output_filename='Precursor_Formula_Ranking_Result_CHONSP.xlsx'; % file name of the formula ranking result
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create a progress bar
hdl = findobj('tag','waitfig');
if ~isempty(hdl) %delete the old progress bar if exist
    delete(hdl);
end
% The output file columns and their variable types:
varNames = {'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference','Rank','Original_Candidates',...
    'Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff',...
        'Peak_Dist_SD','Comp_Time','Note'};
varTypes = {'double','string','double','string','double','double','double','double',...
    'double','double','double','double','double','double',...
        'double','double','string'};
% Find max mass in the peak list
if ~isnumeric(MS2Table.precursor_mz)
    maxmass=max(str2double(MS2Table.precursor_mz)); % max molecular weight in the spectra
else
    maxmass=max(MS2Table.precursor_mz); % max molecular weight in the spectra
end
% Reduce the database according to the max mass
load PubChemMetabolite.mat;
useid=mass <= (maxmass+2.5);
elemnum=uint8(elemnum(useid,:));
formula=formula(useid);
mass=mass(useid);
% convert the input peak list from a table to a struct
spectra=table2struct(MS2Table);
spectnum=size(MS2Table,1); % total number of spectra
% construct a table to record the identification information
spect = table('Size',[spectnum*opt.MaxRankNumber length(varNames)],'VariableTypes',varTypes,'VariableNames',varNames);
% initialize the progress bar
h = waitbar(0, 'Waiting...','tag','waitfig');
% start the identification for each MS/MS spectrum
starttime=cputime;
rec_count=0;
for i=1:spectnum % for each spectrum
    msg=SingleSpectrumFormulaRanking(i,elemnum,formula,mass,spectra(i),opt);
    reclen=size(msg,1);
    if reclen>0
        spect((rec_count+1):(rec_count+reclen),:)=msg;
        rec_count=rec_count+reclen;
    end
    ratio=1.0*i/spectnum;
    waitbar(ratio, h, ['In the given dataset: ',num2str(ratio*100,'%6.2f'),'% (',num2str(i),'/',num2str(spectnum),') finished']);
end
totaltime=cputime-starttime
% write FSM results
writetable(spect,output_filename,'Sheet','formula_ranking_result');
delete(h);
