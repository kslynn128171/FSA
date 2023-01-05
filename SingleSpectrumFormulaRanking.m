function msg=SingleSpectrumFormulaRanking(spectidx, elemnum, formula, mass, spectrum, opt)
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
% varNames = {'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference','Rank','Original_Candidates',...
%     'Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff',...
%         'Peak_Dist_SD','Comp_Time','Note'};
% 
% Author: Ke-Shiuan Lynn Ph.D.
% Assistant Professor
% Department of Mathematics
% Fu-Jen Catholic University
% Email: 128171@mail.fju.edu.tw
% Final Update: Apr. 20, 2022
%--------------------------------------------------------------------------
msg=cell(opt.MaxRankNumber,17);
msg_count=0;
% Extract m/z and abundance values from the "ms2_peak" column
temp = strsplit(spectrum.ms2_peak,';');
tempnum = length(temp);
if rem(tempnum,2) ~= 0 % check the ms2_peak and see if the m/z and abundance come with pairs
    disp(['The format of the ms2_peaks is incorrect in spectrum #',num2str(spectidx),'. The m/z and abundance values should appear in pairs, but ',num2str(tempnum),' values are found.']);
    msg=[];
    return;
end
num=tempnum/2;
terms = reshape(temp,2,num);
mz=str2double(terms(1,:)); % m/z value
ab=str2double(terms(2,:)); % abundance
clear terms
nop_org=length(mz); % number of peaks in the spectrum
% convert to monoisotopic mass based on the ion mode
mode=lower(spectrum.ionization_mode); % the mode of the spectrum
if contains(mode,'os')
    mz=mz-1.007276; % positive mode, perform protonation the peak m/z
    if ~isnumeric(spectrum.precursor_mz)
        precursor=str2double(spectrum.precursor_mz)-1.007276;
    else
        precursor=pectrum.precursor_mz-1.007276;
    end
else
    mz=mz+1.007276; % negative mode, perform deprotonation the peak m/z
    if ~isnumeric(spectrum.precursor_mz)
        precursor=str2double(spectrum.precursor_mz)+1.007276;
    else
        precursor=pectrum.precursor_mz+1.007276;
    end
end
tic; % resume the timer
% ------------------------------------------
% find the mass of fragments in a spectrum
%-------------------------------------------
% compute the SD of inter-peak distances
mzorg=[precursor mz(mz<precursor)];
mzorg=sort(mzorg,'desc');
if length(mzorg) == 1
    Peak_Dist_SD=0;
else
    Peak_Dist_SD=std(mzorg(1:end-1)-mzorg(2:end));
end
% remove peaks with higher mass than the precursor
lobd=precursor-2;
remidx=mz >= lobd; % peak indices to be removed (they are not fragments)
% keep only the precursor and the fragments in the spectrum
if ~any(remidx) % all peaks are smaller than the lower bound (lobd)
    mz=[precursor mz];
    ab=[max(ab) ab];
else % There exist peaks whose m/z values are larger than the lower bound (lobd)
    if any(~remidx) % there are fragments in the spectrum
        mz=[precursor mz(~remidx)]; % m/z values after irrelavant peak removal
        ab=[max(ab(~remidx)) ab(~remidx)]; % abundance values after irrelevant peak removal
    else % no fragment in the spectrum
        massdiff=abs(mass-precursor);
        [sorteddiff,sidx]=sort(massdiff);
        if opt.MatchTolerance < 1
            delta=opt.MatchTolerance;
        else
            delta=opt.MatchTolerance*precursor/1e6;
        end
        is_qualified=sorteddiff<=delta;
        if any(is_qualified) % there exists formula candidates
            comp_time=toc;
            if sum(is_qualified) == 1 % Only one formula candidate is found.
                score=max(0,(delta-sorteddiff(1))/delta);
                msg_count=msg_count+1;
                note='Single candidate exists in the given tolerance.';
                % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
                % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
                msg(msg_count,:)={spectidx,mode,precursor,formula{sidx(1)},mass(sidx(1)),sorteddiff(1),...
                    1,1,nop_org,0,0,score,sorteddiff(1),-1,Peak_Dist_SD,comp_time,note};
            else % there exists multiple formula candidates
                note='No fragment is found. Conventional mass matching is performed.';
                Top2_Score_Diff=(sorteddiff(2)-sorteddiff(1))/delta;
                for i=1:min(sum(is_qualified),opt.MaxRankNumber)
                    msg_count=msg_count+1;
                    % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
                    % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
                    score=max(0,(delta-sorteddiff(i))/delta);
                    msg(msg_count,:)={spectidx,mode,precursor,formula{sidx(i)},mass(sidx(i)),sorteddiff(i),...
                        i,sum(is_qualified),nop_org,0,0,score,sorteddiff(i),Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
                end
            end
        else % No formula candidate can be found under the given tolerance
            note='No formula can be found under the given tolerance.';
            msg_count=msg_count+1;
            % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
            % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
            msg(msg_count,:)={spectidx,mode,precursor,'none',-1,-1,...
                -1,0,nop_org,0,0,-1,-1,-1,Peak_Dist_SD,toc,note};
        end
        msg=msg(1:msg_count,:);
        return
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
% -----------
% start FSA
% -----------
% find possible formula for each m/z values in the spectrum
idrec=cell(nop,1);
difmtx=cell(nop,1);
for j=1:nop
    % adjust match tolerance according to the specified unit
    if opt.MatchTolerance < 1
        delta=opt.MatchTolerance;
    else
        delta=opt.MatchTolerance*mz(j)/1e6;
    end
    upbd=mz(j)+delta;
    lobd=mz(j)-delta;
    tf=(mass <= upbd) & (mass >= lobd);
    idx=find(tf);
    if j == 1 % check for the precursor
        uid=all(elemnum(idx,1:2),2); % must contain both C and H in precursor
        idx=idx(uid);
        if isempty(idx) % No formula is found for the precursor 
            note='No formula can be found containing both C and H.';
            msg_count=msg_count+1;
            % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
            % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
            msg(msg_count,:)={spectidx,mode,precursor,'none',-1,-1,...
                -1,0,nop_org,0,0,-1,-1,-1,Peak_Dist_SD,toc,note};
            msg=msg(1:msg_count,:);
            return;
        end
        if length(idx)==1
            massdiff=abs(mass(idx)-precursor);
            msg_count=msg_count+1;
            score=max(0,(delta-massdiff)/delta);
            note='Single candidate exists in the given tolerance.';
            % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
            % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
            msg(msg_count,:)={spectidx,mode,precursor,formula{idx},mass(idx),massdiff,...
                1,1,nop_org,nop,-1,score,massdiff,-1,Peak_Dist_SD,toc,note};
            msg=msg(1:msg_count,:);
            return;
        end
    end
    idrec{j}=idx; % record fomula indices for each peak
    difmtx{j}=abs(mass(idx)-mz(j)); % record the mass differences of the candidates of the j-th peak
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
    comp_time=toc;
    note='No mother-daughter relation is found for all the precursor formula. Conventional mass matching is performed.';
    idx=idrec{1}; % indices of the precursor candidate
    [sorteddiff,sidx]=sort(abs(difmtx{1})); % sort the mass differences of the precursor candidates
    if min(length(idx),opt.MaxRankNumber) > 1
        Top2_Score_Diff=(sorteddiff(2)-sorteddiff(1))/delta;
    else
        Top2_Score_Diff=-1;
    end
    for i=1:min(length(idx),opt.MaxRankNumber)
        msg_count=msg_count+1;
        % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
        % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
        tempscore=max(0,(delta-sorteddiff(i))/delta);
        msg(msg_count,:)={spectidx,mode,precursor,formula{idx(sidx(i))},mass(idx(sidx(i))),sorteddiff(i),...
            i,length(idx),nop_org,nop,nof,tempscore,sorteddiff(i),Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    end
    msg=msg(1:msg_count,:);
    return
end
totalscore=nof*(nof-1)/2; % the highest score that could possiblly be
difsum=zeros(noc,1);
massid=zeros(nop,noc);
mincnum=zeros(noc,1);
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
        if (cnum(k)==1) 
            massid(k,j)=find(uid);
        end
    end
    if ~isempty(cnum(cnum>0))
        mincnum(j)=min(cnum(cnum>0));
    end
    uidrec=uidrec(2:end); % remove the first record (unused) in the formula indices
    cnum=cnum(2:end); % remove the first record (unused) in the candidate number
    massuid=massid(2:end,j); % remove the first record (unused) in the mass indices
    uidrec=cell2mat(uidrec(cnum==1)); % keep only characteristic fragments
    massuid=massuid(cnum==1);
    nor_cur=length(uidrec); % number of fragments whose formula can be found
    % compute the match score of each combination
    [k,l]=ind2sub([nor_cur,nor_cur],nonzeros(triu(reshape(1:(nor_cur*nor_cur), [nor_cur,nor_cur]),1)));
    tempscore=sum(all(elemnum(uidrec(k),:) >= elemnum(uidrec(l),:),2));
    score(j)=(tempscore+nor_cur)/totalscore;
    precrec(j)=idrec{1}(j);
    nor(j)=nor_cur;
    % compute the mass differences of this combination
    tempmtx=difmtx(2:end); % keep the mass differences of the fragments
    tempmtx=tempmtx(cnum==1); % keep the mass differences of the characteristic fragments
    difvec=cellfun(@(x,y)x(y),tempmtx,num2cell(massuid)); % extract the mass differences of the characteristic fragments
    difsum(j)=sum(difvec)+difmtx{1}(j); % add the mass difference of the precursor
end
comp_time=toc; % total computation time
if all(score==0) % No mother-daughter relation is found for all the precursor formula
    note='No mother-daughter relation is found for all the precursor formula. Conventional mass matching is performed.';
    idx=idrec{1};
    [sorteddiff,sidx]=sort(difmtx{1}); %mass differences of the precursor formulas
    if min(length(idx),opt.MaxRankNumber) > 1
        Top2_Score_Diff=(sorteddiff(2)-sorteddiff(1))/delta;
    else
        Top2_Score_Diff=-1;
    end
    for i=1:min(length(idx),opt.MaxRankNumber)
        msg_count=msg_count+1;
        % 'No','Ionization_Mode','Precursor_Mass','Formula','Formula_Mass','Mass_Difference',...
        % 'Rank','Original_Candidates','Peak_Num','Frag_Num','Frag_w_Formula','Score','Difference_Sum','Top2_Sore_Diff','Peak_Dist_SD','Comp_Time','Note';
        tempscore=max(0,(delta-sorteddiff(i))/delta);
        msg(msg_count,:)={spectidx,mode,precursor,formula{idx(sidx(i))},mass(idx(sidx(i))),sorteddiff(i),...
            i,length(idx),nop_org,nop,nof,tempscore,sorteddiff(i),Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
    end
else
    idx=idrec{1};
    [ss,sidx]=sort(score,'descend'); % sort the score in a descending order
    Top2_Score_Diff=ss(1)-ss(2); % the difference of the top2 scores
    [uscore,~,uid]=unique(-score); % find unique scores
    uscore=-uscore; % scores descendent sorting
    % find out the rank of the correct formula after adjusting for mass difference 
    rank=1;
    % disp the top-opt.MaxRankNumber precursor formulas
    note='successful ranked';
    for i=1:length(uscore) 
        tienum=sum(uid==i);
        if tienum==1 % no candidate has the same score
            msg_count=msg_count+1;
            msg(msg_count,:)={spectidx,mode,precursor,formula{idx(sidx(i))},mass(idx(sidx(i))),difmtx{1}(sidx(i)),...
                rank,length(idx),nop_org,nop,nor(sidx(i)),uscore(i)+1,-1,Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
            rank=rank+1;
            if rank > opt.MaxRankNumber
                break;
            end
        else % There are ties in the candidates
            [sdif,didx]=sort(difsum(uid==i)); % sort the mass difference
            tieidx=find(uid==i);
            for k=1:length(didx)
                useid=sidx(tieidx(didx(k)));
                msg_count=msg_count+1;
                msg(msg_count,:)={spectidx,mode,precursor,formula{idx(useid)},mass(idx(useid)),difmtx{1}(useid),...
                    rank,length(idx),nop_org,nop,nor(useid),score(useid)+1,...
                    sdif(k),Top2_Score_Diff,Peak_Dist_SD,comp_time,note};
                rank=rank+1;
                if rank > opt.MaxRankNumber
                    break;
                end
            end
            if rank > opt.MaxRankNumber
                break;
            end
        end
    end
end
msg=msg(1:msg_count,:);