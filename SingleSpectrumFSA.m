% -------------------------------------------------
% Single MS/MS spectrum Identification using FSA
% -------------------------------------------------
% Author: Ke-Shiuan Lynn Ph.D.
% Email: 128171@mail.fju.edu.tw
% Final Update: Nov. 12, 2021
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
Match_Tol=0.0012; %recommanded values: QTOF: 0.0012, ITFT: 0.0011, QFT: 0.0008
MaxQualifiedResult=5;
%---------------------------------------------------------------------------
% load PubChem database
load PubChemDatabase
% initialize parameters
TINY=0.0;
note='successful ranked';
Top2_Score_Diff=-1;
comp_time=-1;
is_finished=false;
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
    disp('No mother-daughter relation is found for all the precursor formula.');
    disp('--------------------- Identification Result (Mass Match)-----------------')
    idx=idrec{1};
    massdiff=mass(idx)-precursor;
    [sorteddiff,sidx]=sort(massdiff);
    for i=1:length(idx)
        fprintf('Formula: %15s, Mass: %10.4f, Difference: %10.5e\n',formula{idx(sidx(i))},mass(idx(sidx(i))),sorteddiff(i));
    end
else
    idx=idrec{1};
    [ss,sidx]=sort(score,'descend'); % sort the score in a descending order
    Top2_Score_Diff=ss(1)-ss(2); % the difference of the top2 scores
    uscore=sort(unique(score),'desc'); % score sorting
    bestscore=max(score); % the best score among all formula cansidates
    CharFragPerc=nor/nof; % percentage of characteristic fragment to total fragments for all formula candidates
    MaxCF=max(CharFragPerc); % max. percentage of characteristic fragment
    ranknum=find(uscore==score(1)); % score rank (same scores are ignored)
    ranking=sum(score>=score(1)); % score rank (same scores are accounted for)
    is_top1_uniq=sum(score==max(score(1)))==1;
    nor_rec=nor(1);
    disp('--------------------- Identification Result (FSA)-----------------')
    for i=1:length(score)
        fprintf('Formula: %15s, Mass: %10.4f, Score: %10.4f\n',formula{idx(sidx(i))},mass(idx(sidx(i))),ss(i));
    end
    disp(['SD of neighboring peak distance = ',num2str(Peak_Dist_SD,'%10.4f')]);
    disp(['Characteristic fragment ratio = ',num2str(CharFragPerc(sidx(1)),'%10.4f')]);
    disp(['Top-2 score difference = ',num2str(Top2_Score_Diff,'%10.4f')]);
end
