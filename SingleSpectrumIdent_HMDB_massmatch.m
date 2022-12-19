function msg=SingleSpectrumIdent_HMDB_massmatch(spectidx, elemnum, formula, mass, spectrum, opt)
% initial parameters
note='successful ranked';
Top2_Score_Diff=-1;
comp_time=-1;
is_reduced=false;
instrument_type=strrep(spectrum.instrument_type,',',';');
% split m/z and abundance values in a spectral file
temp = strsplit(spectrum.ms2_peak,{'*',';'});
num = length(temp);
if rem(num,4)==0
    terms = reshape(temp,4,num/4);
else
    terms = reshape(temp(1:(floor(num/4)*4)),4,floor(num/4));
end
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
if ~isnumeric(spectrum.monisotopic_molecular_weight)
    precursor=str2double(spectrum.monisotopic_molecular_weight); % precursor mass
else
    precursor=spectrum.monisotopic_molecular_weight; % precursor mass
end
mzdiff=abs(mz-precursor);
[mindiff,minidx]=min(mzdiff);
% compute the SD of inter-peak distances
mzorg=[precursor mz(mz<precursor)];
mzorg=sort(mzorg,'desc');
if length(mzorg) == 1
    Peak_Dist_SD=0;
else
    Peak_Dist_SD=std(mzorg(1:end-1)-mzorg(2:end));
end
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
        if mindiff < opt.MatchTolerance
            mz=mz(minidx);
        else
            mz=precursor;
        end
        ab=1.0;
    end
end
% normalize abundance 
ab=ab/max(ab)*100;
% find the top-insensed peaks
sab=sort(ab,'desc');
if opt.MaxPeak > 1
    rid=(ab>=sab(min(length(ab),opt.MaxPeak))); % a fixed peak number
else
    rid=(ab>=sab(min(length(ab),floor(opt.MaxPeak*precursor)))); % a ratio of the precurcor mass
end
ab=ab(rid);mz=mz(rid);
% remove tiny peaks
kid=ab > opt.TINY;
mz=mz(kid);
% sort the remaining peaks in descending m/z order
mz=sort(mz,'desc');
nop=length(mz); % number of peaks in the spectrum
% find possible formula for each m/z values in the spectrum
idrec=cell(nop,1);
difmtx=cell(nop,1);
for j=1:nop
    % adjust match tolerance according to the specified unit
    if opt.MatchTolerance < 1 % use Da as unit of the tolerance
        delta=opt.MatchTolerance;
    else % use PPM as unit of the tolerance
        if mz(j) <= opt.mass
            delta=opt.tol;
        else % relax the tolerance for small mass
            delta=opt.MatchTolerance*mz(j)/1e6;
        end
    end
    upbd=mz(j)+delta;
    lobd=mz(j)-delta;
    tf=(mass <= upbd) & (mass >= lobd);
    idx=find(tf);
    idrec{j}=idx; % record fomula indices for each peak
    difmtx{j}=mass(idx)-mz(j); % record the mass differences of the candidates
    if j == 1 % check for the precursor
        uid=all(elemnum(idx,1:2),2); % must contain both C and H in precursor
        idx=idx(uid);
        if isempty(idx) % No formula is found for the precursor 
            note='No formula is found for the precursor.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                -1,-1,-1,nop_org,nop,-1,-1,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
        tf=strcmp(formula(idx),spectrum.chemical_formula);
        idx=[idx(tf);idx(~tf)]; % move the correct answer to the top
        idrec{j}=idx; % record fomula indices for each peak
        difmtx{j}=mass(idx)-mz(j); % record the mass differences of the candidates
        score=max(0,(opt.MatchTolerance-min(difmtx{j}))/opt.MatchTolerance);
        if ~any(tf) % No precursor formula matches with the answer
            note='No precursor formula matches with the answer.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                -1,-1,length(idx),nop_org,nop,-1,score,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
        if length(idx)==1
            note='Single candidate matches with the answer.';
            msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
                1,1,1,nop_org,nop,-1,score,-1,-1,-1,-1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            return;
        end
    end 
end
%-----------------------------------------------------------
% find mother-daughter relations for each precursor formula
%-----------------------------------------------------------
noc=length(idrec{1}); % number of formula candidates for the precursor
nor=zeros(noc,1); % number of fragments with mother-daughter relationships (for a certain combination)
score=zeros(noc,1); % percentage of mother-daughtor relationships can be found (for a certain combination)
nof=sum(~cellfun(@isempty,idrec)); % number of fragment whose formula can be found
if nof == 1
    note='No mother-daughter relation is found for all the precursor formula.';
    [smass,sidx]=sort(abs(difmtx{1}));
    [usmass,~]=unique(smass);
    % use the score from mass match
    bestscore=max(0,(opt.MatchTolerance-smass(1))/opt.MatchTolerance);
    is_top1_uniq=smass(1)==smass(2);
    ranking=find(sidx==1);
    rank=find(usmass==smass(ranking),1,'first');
    msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
        rank,ranking,noc,nop_org,nop,nor(1),bestscore,-1,is_top1_uniq,-1,-1,-1,...
        Top2_Score_Diff,Peak_Dist_SD,0,note};
    return
end
totalscore=nof*(nof-1)/2; % the highest score that could possiblly be
difsum=zeros(noc,1);
massid=zeros(nop,noc);
tic; % resume the timer
for j=1:noc % for each precursor formula candidate
    vecp=elemnum(idrec{1}(j),:); % element composition in the precursor formula
    cnum=ones(nop,1);
    uidrec=cell(nop,1);
    for k=2:nop % make sure that all the fragment components are a subset of those of the precursor
        idx=idrec{k}; % indices of formulas for the kth fragment
        matf=elemnum(idx,:); % fragment matrix
        matp=repmat(vecp,length(idx),1); % parent matrix
        tf=matp >= matf; % logical matrix indicating whether a subset is found
        uid=all(tf,2); % logical value showing which formulas that possess MDR with precursor
        cnum(k)=sum(uid); % record number of satisfied candidates
        uidrec{k}=(idx(uid))'; % record the indices of the satisfied candidates
        if (cnum(k)==1) 
            massid(k,j)=find(uid);
        end
    end
    combnum=prod(cnum(cnum>0)); % total number of combination paths to be examined
    if (j == 1) &&  (combnum > 1) % record this info if the formula is the correct answer
        is_reduced=true; % the combination of the correct formula contains peak(s) with multiple formulas
    end
    uidrec=uidrec(2:end); % remove the first record (unused) in the formula indices
    cnum=cnum(2:end); % remove the first record (unused) in the candidate number
    massuid=massid(2:end,j); % remove the first record (unused) in the mass indices
    uidrec=cell2mat(uidrec(cnum==1)); % keep only characteristic fragments
    massuid=massuid(cnum==1);
    nor_cur=length(uidrec); % number of characteristic fragments
    % compute the match score of each combination
    [k,l]=ind2sub([nor_cur,nor_cur],nonzeros(triu(reshape(1:(nor_cur*nor_cur), [nor_cur,nor_cur]),1)));
    tempscore=sum(all(elemnum(uidrec(k),:) >= elemnum(uidrec(l),:),2)); % number of MDRs among the fragments
    score(j)=(tempscore+nor_cur)/totalscore; % add the MDRs with the precursor and normalize the score
    nor(j)=nor_cur;
    % compute the mass differences of this combination
    tempmtx=difmtx(2:end); % keep the mass differences of the fragments
    tempmtx=tempmtx(cnum==1); % keep the mass differences of the characteristic fragments
    difvec=cellfun(@(x,y)x(y),tempmtx,num2cell(massuid)); % extract the mass differences of the characteristic fragments
    difsum(j)=(sum(abs(difvec))+abs(difmtx{1}(j)))/(length(difvec)+1); % add the mass difference of the precursor
end
comp_time=toc; % total computation time
if all(score==0) % No mother-daughter relation is found for all the precursor formula
    note='No mother-daughter relation is found for all the precursor formula.';
    [smass,sidx]=sort(abs(difmtx{1}));
    [usmass,~]=unique(smass);
    % report the score from mass match
    bestscore=max(0,(opt.MatchTolerance-smass(1))/opt.MatchTolerance);
    is_top1_uniq=smass(1)==smass(2);
    ranking=find(sidx==1);
    rank=find(usmass==smass(ranking),1,'first');
    msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
        rank,ranking,noc,nop_org,nop,nor(1),bestscore,-1,is_top1_uniq,-1,-1,-1,...
        Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    return
% elseif score(1) == 0 % No mother-daughter relation is found for the correct answer
%     note='No mother-daughter relation is found for the precursor.';
%     msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
%         0,0,noc,nop_org,nop,nor(1),0,is_reduced,0,0,0,0,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
%     return
end
ss=sort(score,'descend'); % sort the score in a descending order
Top2_Score_Diff=ss(1)-ss(2); % the difference of the top2 scores
[uscore,~,uid]=unique(-score);
uscore=-uscore;
bestscore=max(score)+1; % the best score among all formula cansidates
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
ranking=sum(score>=score(1)); % score rank (same scores are accounted for)
% find out the rank of the correct formula after adjusting for mass difference 
rank=0;
for i=1:length(uscore)
    tienum=sum(uid==i);
    if uscore(i)==score(1) % the candidates with the same score as the correct precursor formula
        if tienum==1 % no candidate has the same score
            rank=rank+1;
            break;
        else
            sdif=sort(difsum(uid==i)); % sort the mass difference
            rid=find(sdif==difsum(1)); % find the indices of the sorted differences that matches that of the correct formula
            rank=rank+rid(1); % use the rank of the mass difference to break the ties
            break;
        end
    else % candidates that have scores better than the correct precursor formula
        rank=rank+tienum;
    end
end
is_top1_uniq=sum(score==max(score(1)))==1;
nor_rec=nor(1);
msg={spectidx,spectrum.name,instrument_type,mode,precursor,spectrum.chemical_formula,...
    rank,ranking,noc,nop_org,nop,nor_rec,bestscore,is_reduced,is_top1_uniq,CharFragPerc(1),isMAXCF,isMAXCF_uniq,...
    Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
