% 96330261 pubmed records -> 94215477 CHONSP+FClBrI records -> 2613203
% unique records of which 1970 records have their masses corrected
clearvars
clc
close all
PBTable=readtable('PubChem_metabolites.csv');
formula=table2cell(PBTable(:,1));
mass=table2array(PBTable(:,2));
clear PBTable
recnum=size(formula,1);
KeepElem={'C','H','O','N','S','P'};%,'F','Cl','Br','I'};
elemnum=zeros(recnum,length(KeepElem));
SelectIdx=false(recnum,1);
count=0;
for i=1:recnum
    str=formula{i};
    [EleList,Trash,EleEnd]=regexp(str,['[','A':'Z','][','a':'z',']?'],'match');
    [Num,NumStart]=regexp(str,'\d+','match');
    NumList=ones(size(EleList));
    Index=ismember(EleEnd+1,NumStart);
    NumList(Index)=cellfun(@str2num,Num);
    target = arrayfun(@(x)char(EleList(x)),1:numel(EleList),'uni',false);
    index = cellfun(@(a) strcmpi(a,KeepElem),target,'uniform',false);
    if all(cellfun(@any,index))
        tempmat=cell2mat(index');
        tempvec=NumList*tempmat;
        count=count+1;
        elemnum(count,:)=tempvec;
        SelectIdx(i)=true;
    end
    if rem(i,1000000)==0
        disp(i)
    end
end
% remove the additional space
elemnum=elemnum(1:count,:);
formula=formula(SelectIdx);
mass=mass(SelectIdx);
% remove repeated record
[formula,uid]=unique(formula);
elemnum=elemnum(uid,:);
mass=mass(uid);
% remove element numbers that are greater than 255
revid=find(elemnum>255);
[row,col]=ind2sub(size(elemnum),revid);
elemnum(row,:)=[];
formula(row)=[];
mass(row)=[];
%% correct mass
element_mass=[12.0000000000;1.0078250321;15.9949146221;14.0030740052;31.9720706900;...
    30.9737615100];%18.9984032000;34.9688527100;78.9183376000;126.9044680000];
massnew=double(elemnum)*element_mass;
diff=abs(massnew-mass);
id=find(diff>=0.0005);
length(id)
mass=massnew;
elemnum=uint8(elemnum);
save PubChemMetabolite.mat mass formula elemnum