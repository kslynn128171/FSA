% -----------------------------------------------------------------
% Precursor formula ranking using FSA for a single MS/MS peak list
% -----------------------------------------------------------------
% Author: Ke-Shiuan Lynn Ph.D.
% Email: 128171@mail.fju.edu.tw
% Final Update: Apr. 20, 2022
%--------------------------------------------------------------------------
% clear memory
clearvars
clc
close all
% -------------------------------------------------------------------------
% Please fill in the following required information
precursor_mz=245.2341;
mz=[245.2341, 171.1502, 129.1395, 112.1128, 100.0760, 84.0812];
ab=[56.8106, 58.1811, 65.6561, 90.0748, 100.0000, 51.4535];
mode='positive';
Match_Tol=0.001; 
MaxQualifiedResult=5;
%---------------------------------------------------------------------------
% load PubChem database
load PubChemMetabolite
% initialize parameters
TINY=0.0; % threshold to remove peaks with small abundance
% ------------------------------------------------------------------------
% Preprocessing the m/z values
nop_org=length(mz); % number of peaks in the spectrum
% adjust mode for the m/z values
if contains(lower(mode),'os')
    mz=mz-1.007276; % positive mode, perform protonation the peak m/z
    precursor=precursor_mz-1.007276;
else
    mz=mz+1.007276; % negative mode, perform deprotonation the peak m/z
    precursor=precursor_mz+1.007276;
end
% ------------------------------------------
% find the mass of fragments in a spectrum
%-------------------------------------------
mzorg=mz;
% compute the SD of inter-peak distances
mztmp=[precursor mz(mz<precursor)];
mztmp=sort(mztmp,'desc');
Peak_Dist_SD=std(mztmp(1:end-1)-mztmp(2:end));
% remove peaks with higher mass than the precursor
lobd=precursor-2;
remidx=mz >= lobd; % peak indices to be removed (they are not fragments)
% keep only the precursor and the fragments in the spectrum
if ~any(remidx) % all peaks are smaller than the lower bound (lobd)
    mz=[precursor mz];
    ab=[max(ab) ab];
else
    if any(~remidx) % there are fragments in the spectrum
        mz=[precursor mz(~remidx)]; % m/z values after irrelavant peak removal
        ab=[max(ab(~remidx)) ab(~remidx)]; % abundance values after irrelevant peak removal
    else % no fragment in the spectrum
        massdiff=abs(mass-precursor);
        [sorteddiff,sortidx]=sort(massdiff);
        disp('No fragment is found in the spectrum.');
        disp(['    Number of peaks: ',num2str(nop_org);]);
        disp('    Number of fragments: 0');
        disp('--------------------- Identification Result (Mass Match)-----------------')
        is_qualified=sorteddiff<=Match_Tol;
        if any(is_qualified)
            for i=1:sum(is_qualified)
                fprintf('Formula: %15s, Mass: %10.4f, Difference: %10.5e\n',formula{sortidx(i)},mass(i),sorteddiff(i));
            end
        else
            disp(['No qualified formula can be found: minimum difference ',num2str(sorteddiff(1),'%10.5e'),' required tolerance: ',num2str(Match_Tol,'10.5e')]);
        end
        return
    end
end
% normalize abundance 
ab=ab/max(ab)*100;
% remove tiny peaks
kid=ab > TINY;
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
    % adjust match tolerance if the molecular mass is small
    upbd=mz(j)+Match_Tol;
    lobd=mz(j)-Match_Tol;
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
            disp('--------------------- Identification Result (FSA)-----------------------')
            disp('No formula can be found for the precursor.');
            return
        end
        if length(idx)==1
            massdiff=abs(mass(idx)-precursor);
            disp('--------------------- Identification Result (FSA)-----------------------')
            disp('Single candidate matches with the answer:')
            fprintf('Formula: %15s, Mass: %10.4f, Difference: %10.5e\n',formula{idx},mass(idx),massdiff);
            return
        end
    end
    idrec{j}=idx; % record fomula indices for each peak
    difmtx{j}=abs(mass(idx)-mz(j)); % record the mass differences of the candidates
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
    disp('No mother-daughter relation is found for all the precursor formula! Mass match is performed');
    disp('--------------------- Identification Result (Mass Match)-----------------')
    idx=idrec{1};
    massdiff=mass(idx)-precursor;
    [sorteddiff,sidx]=sort(massdiff);
    for i=1:length(idx)
        fprintf('Formula: %15s, Mass: %10.4f, Difference: %10.5e\n',formula{idx(sidx(i))},mass(idx(sidx(i))),sorteddiff(i));
    end
    return
end
totalscore=nof*(nof-1)/2; % the highest score that could possiblly be
difsum=zeros(noc,1);
massid=zeros(nop,noc);
mincnum=zeros(noc,1);
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
disp(['Computation Time: ',num2str(comp_time),' sec']);
if all(score==0) % No mother-daughter relation is found for all the precursor formula
    disp('No mother-daughter relation is found for all the precursor formula.');
    disp('--------------------- Identification Result (Mass Match)-----------------')
    idx=idrec{1};
    [sorteddiff,sidx]=sort(difmtx{1}); %mass differences of the precursor formulas
    for i=1:min(MaxQualifiedResult,length(score))
        fprintf('Formula: %15s, Mass: %10.4f, Difference: %10.5e\n',formula{idx(sidx(i))},mass(idx(sidx(i))),sorteddiff(i));
    end
else
    idx=idrec{1};
    [ss,sidx]=sort(score,'descend'); % sort the score in a descending order
    Top2_Score_Diff=ss(1)-ss(2); % the difference of the top2 scores
    [uscore,uid]=sort(unique(score),'desc'); % score sorting
    bestscore=max(score); % the best score among all formula cansidates
    %ranknum=find(uscore==score(1)); % score rank (same scores are ignored)
    ranking=sum(score>=score(1)); % score rank (same scores are accounted for)
    % find out the rank of the correct formula after adjusting for mass difference 
    disp('--------------------- Identification Result (FSA)-----------------')
    disp(['SD of neighboring peak distance = ',num2str(Peak_Dist_SD,'%10.4f')]);
    disp(['Top-2 score difference = ',num2str(Top2_Score_Diff,'%10.4f')]);
    rank=1;
    % disp the top-MaxQualifiedResult precursor formulas
    for i=1:length(uscore) 
        tienum=sum(uid==i);
        if tienum==1 % no candidate has the same score
            fprintf('Rank#%d, Formula: %15s, Mass: %10.4f, Score: %10.4f\n',rank,formula{idx(sidx(rank))},mass(idx(sidx(rank))),score(sidx(rank)));
            disp(['Characteristic fragment ratio = ',num2str(nor(sidx(rank))/nof,'%10.4f')]);
            rank=rank+1;
            if rank > MaxQualifiedResult
                break;
            end
        else
            [sdif,didx]=sort(difsum(uid==i)); % sort the mass difference
            tieidx=find(uid==idrec);
            tempformula=formula(idx(sidx(tieidx)));
            tempmass=mass(idx(sidx(tieidx)));
            for k=1:length(didx)
                fprintf('Rank#%d, Formula: %15s, Mass: %10.4f, Score: %10.4f, Sum of Mass Differences: %10.4f\n',rank,tempformula{didx(k)},tempmass(didx(k)),score(sidx(rank)),sdif(k));
                disp(['Characteristic fragment ratio = ',num2str(nor(sidx(rank))/nof,'%10.4f')]);
                rank=rank+1;
                if rank > MaxQualifiedResult
                    break;
                end
            end
            if rank > MaxQualifiedResult
                break;
            end
        end
    end
end