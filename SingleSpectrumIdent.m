function msg=SingleSpectrumIdent(spectidx, elemnum, formula, mass, spectrum, opt)
% SingleSpectrumIdent Rank the precursor formula candidates using formula
% subset analysis.
% input parameters: 
% spectidx: spectrum index
% elemnum: the matrix that contains C, H, O, N, S, and P element numbers for compounds in PubChem
% formula: the cell array that contains compound formulas in PubChem
% mass: the vector that contains monoisotope mass for compounds in PubChem
% spectrum: the MS/MS spectrum to be identified
% opt: parameter options
% The output file columns are: 
% varNames = {'No','Name','Instrument_Type','Ionization_Mode','Mass','Formula','Rank_Reduced','Original_Rank','Original_Candidates','Peak_Num',...
%         'Frag_Num','Frag_w_Formula','Score','Is_Reduced','Is_Top1_Uniq','Char_Frag_Percent','Is_CFP_Max','Is_CFP_Max_Uniq','Top2_Sore_Diff',...
%         'Peak_Dist_SD','Comp_Time','Note'};
% 
% Author: Ke-Shiuan Lynn Ph.D.
% Assistant Professor
% Department of Mathematics
% Fu-Jen Catholic University
% Email: 128171@mail.fju.edu.tw
% Final Update: Nov. 12, 2021
%--------------------------------------------------------------------------

% initial parameters
note='successful ranked';
Top2_Score_Diff=-1;
comp_time=-1;
is_reduced=false;
instrument_type=strrep(spectrum.instrument_type,',',';');
% split m/z and abundance values in a spectral file
temp = strsplit(spectrum.ms2_peak,{'*',';'});
num = length(temp);
terms = reshape(temp,4,num/4);
mz=str2double(terms(3,:)); % m/z value
ab=str2double(terms(2,:)); % abundance
clear terms
nop_org=length(mz); % number of peaks in the spectrum
% adjust mode for the m/z values
mode=lower(spectrum.ionization_mode); % the mode of the spectrum
if ~contains(mode,'os')
    mz=mz+1.007276; % negative mode, perform deprotonation the peak m/z
else
    mz=mz-1.007276; % positive mode, perform protonation the peak m/z
end
% ------------------------------------------
% find the mass of fragments in a spectrum
%-------------------------------------------
precursor=str2double(spectrum.monisotopic_molecular_weight); % precursor mass
mzdiff=abs(mz-precursor);
[mindiff,minidx]=min(mzdiff);
% compute the SD of inter-peak distances
mzorg=[precursor mz(mz<precursor)];
mzorg=sort(mzorg,'desc');
Peak_Dist_SD=std(mzorg(1:end-1)-mzorg(2:end));
lobd=precursor-2;
cidx=mz >= lobd; % peak indices to be removed (they are not fragments)
ridx=~cidx; % peak indices to be kept (they are most likely fragments)
% keep only the precursor and the fragments in the spectrum
if ~any(cidx) % all peaks are smaller than the lower bound (lobd)
    mz=[precursor mz];
    ab=[max(ab) ab];
else
    if any(ridx) % there are fragments in the spectrum
        mz=[precursor mz(ridx)]; % m/z values after irrelavant peak removal
        ab=[max(ab(ridx)) ab(ridx)]; % abundance values after irrelevant peak removal
    else % no fragment in the spectrum
        if mindiff < 0.01
            mz=mz(minidx);
            ab=1.0;
        else
            note='No close formula to the theoretical precursor.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                -1,-1,-1,nop_org,0,-1,-1,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
    end
end
% normalize abundance 
ab=ab/max(ab)*100;
% remove tiny peaks
kid=ab > opt.TINY;
mz=mz(kid);
% sort the remaining peaks in descending m/z order
mz=sort(mz,'desc');
nop=length(mz); % number of peaks in the spectrum
% find possible formula for each m/z values in the spectrum
idrec=cell(nop,1);
for j=1:nop
    % adjust match tolerance if the molecular mass is small
    upbd=mz(j)+opt.PPM;
    lobd=mz(j)-opt.PPM;
    tf=(mass <= upbd) & (mass >= lobd);
    idx=find(tf);
    if j == 1 % check for the precursor
        uid=false(size(idx));
        for k=1:length(idx) % must contain both C and H in precursor
            if elemnum(idx(k),1)>0 && elemnum(idx(k),2)>0
                uid(k)=true;
            end
        end
        idx=idx(uid);
        if isempty(idx) % No formula is found for the precursor 
            note='No formula is found for the precursor.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                -1,-1,-1,nop_org,nop,-1,-1,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
        tf=strcmp(formula(idx),spectrum.chemical_formula);
        if ~any(tf) % No precursor formula matches with the answer
            note='No precursor formula matches with the answer.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                -1,-1,length(idx),nop_org,nop,-1,-1,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
        idx=[idx(tf);idx(~tf)]; % move the correct answer to the top
        if length(idx)==1
            note='Single candidate matches with the answer.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                1,1,1,nop_org,nop,-1,-1,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
    end
    idrec{j}=idx; % record fomula indices for each peak
end
%-----------------------------------------------------------
% find mother-daughter relations for each precursor formula
%-----------------------------------------------------------
noc=length(idrec{1}); % number of formula candidates for the precursor
nor=zeros(noc,1); % number of fragments with mother-daughter relationships (for a certain combination)
score=zeros(noc,1); % percentage of mother-daughtor relationships can be found (for a certain combination)
precrec=zeros(noc,1); % the precursor formula (for a certain combination)
nof=sum(~cellfun(@isempty,idrec)); % number of fragment whose formula can be found
if nof == 1
    note='No mother-daughter relation is found for all the precursor formula.';
    msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
        -1,-1,noc,nop_org,nop,nor(1),-1,0,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    return
end
totalscore=nof*(nof-1)/2; % the highest score that could possiblly be
tic; % resume the timer
for j=1:noc % for each precursor formula candidate
    vecp=elemnum(idrec{1}(j),:); % element composition in the precursor formula
    cnum=ones(nop,1);
    uidrec=cell(nop,1);
    for k=2:nop % make sure that all the fragment components are a subset of those of the precursor
        idx=idrec{k};
        matf=elemnum(idx,:); % fragment matrix
        matp=repmat(vecp,length(idx),1); % parent matrix
        tf=matp >= matf;
        uid=all(tf,2);
        cnum(k)=sum(uid);
        uidrec{k}=uint32(idx(uid))';
    end
    combnum=prod(cnum(cnum>0)); % total number of combination paths to be examined
    if (j == 1) &&  (combnum > 1) % record this info if the formula is the correct answer
        is_reduced=true; % the combination of the correct formula contains peak(s) with multiple formulas
    end
    uidrec=uidrec(2:nop);
    cnum=cnum(2:nop);
    uidrec=cell2mat(uidrec(cnum==1)); % keep only characteristic fragments
    nor_cur=length(uidrec); % number of fragments whose formula can be found
    % compute the match score of each combination
    [k,l]=ind2sub([nor_cur,nor_cur],nonzeros(triu(reshape(1:(nor_cur*nor_cur), [nor_cur,nor_cur]),1)));
    tempscore=sum(all(elemnum(uidrec(k),:) >= elemnum(uidrec(l),:),2));
    score(j)=(tempscore+nor_cur)/totalscore;
    precrec(j)=idrec{1}(j);
    nor(j)=nor_cur;
end
comp_time=toc; % total computation time
if all(score==0) % No mother-daughter relation is found for all the precursor formula
    note='No mother-daughter relation is found for all the precursor formula.';
    msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
        -1,-1,noc,nop_org,nop,nor(1),-1,is_reduced,0,0,0,0,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    return
elseif score(1) == 0 % No mother-daughter relation is found for the correct answer
    note='No mother-daughter relation is found for the precursor.';
    msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
        0,0,noc,nop_org,nop,nor(1),0,is_reduced,0,0,0,0,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    return
end
ss=sort(score,'descend'); % sort the score in a descending order
Top2_Score_Diff=ss(1)-ss(2); % the difference of the top2 scores
uscore=sort(unique(score),'desc'); % score sorting
bestscore=max(score); % the best score among all formula cansidates
CharFragPerc=nor/nof; % percentage of characteristic fragment to total fragments for all formula candidates
MaxCF=max(CharFragPerc); % max. percentage of characteristic fragment
if CharFragPerc(1) == MaxCF
    isMAXCF = 1; % the percentage of characteristic fragment of the correct formula is the max.
else
    isMAXCF = 0;
end
if sum(CharFragPerc == MaxCF)==1
    isMAXCF_uniq = 1; % the percentage of characteristic fragment of the correct formula is the sole max.
else
    isMAXCF_uniq = 0;
end
ranknum=find(uscore==score(1)); % score rank (same scores are ignored)
ranking=sum(score>=score(1)); % score rank (same scores are accounted for)
is_top1_uniq=sum(score==max(score(1)))==1;
nor_rec=nor(1);
msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
    ranknum,ranking,noc,nop_org,nop,nor_rec,bestscore,is_reduced,is_top1_uniq,CharFragPerc(1),isMAXCF,isMAXCF_uniq,...
    Top2_Score_Diff,Peak_Dist_SD,comp_time,note};

