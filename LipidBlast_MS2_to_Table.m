clearvars
clc
close all
fclose all;
% read positive ms/ms spectrara
fid=fopen('E:\LipidBlast-Full-Release-3\LipidBlast-ASCII-spectra\LipidBlast-neg.msp','r');
spectnum=140000;
name=cell(spectnum,1);
chemical_formula=cell(spectnum,1);
monisotopic_molecular_weight=zeros(spectnum,1);
precursor_MZ=zeros(spectnum,1);
ion=cell(spectnum,1);
ionization_mode=cell(spectnum,1);
peak_num=zeros(spectnum,1);
ms2_peak=cell(spectnum,1);
spectra=table(name,chemical_formula,monisotopic_molecular_weight,precursor_MZ,...
    ion,ionization_mode,peak_num,ms2_peak);
count=1;
emptyline=0;
while ~feof(fid)
    str=fgetl(fid);
    if length(str)<5
        emptyline=emptyline+1;
        if emptyline==1
            count=count+1;
        end
        continue;
    else
        emptyline=0;
    end
    [item,data]=strtok(str,':');
    num=str2double(item);
    switch item
        case 'Name'
            info=strsplit(data(2:end),';');
            spectra.name{count}=info{1};
            spectra.ion{count}=info{2};
            %if info{2}(end)=='-'
                spectra.ionization_mode{count}='negative';
            %else
            %    spectra.ionization_mode{count}='positive';
            %end
        case 'PRECURSORMZ'
            spectra.precursor_MZ(count)=str2double(data(2:end)); % compound precursor m/z
            if strcmpi(spectra.ionization_mode{count},'positive')
                spectra.monisotopic_molecular_weight(count)=spectra.precursor_MZ(count)-1.007276;
            else
                spectra.monisotopic_molecular_weight(count)=spectra.precursor_MZ(count)+1.007276;
            end
        case 'Num Peaks'
            spectra.peak_num(count)=str2double(data(2:end)); % number of peaks in MS2
            tempstr=[];
            for i=1:spectra.peak_num(count)
                str=fgetl(fid);
                info=strsplit(str,' ');
                tempstr=strcat(tempstr,[';0;',num2str(info{2}),';',num2str(info{1}),';0']);
            end
            spectra.ms2_peak{count}=tempstr(2:end);
        case 'Comment'
            info=strsplit(strtrim(data(2:end)),';');
            if length(info{end})<3
                spectra.chemical_formula{count}=info{end-1}; % compound formula
            else
                spectra.chemical_formula{count}=info{end}; % compound formula
            end
    end
end
spectra=spectra(1:count,:);
writetable(spectra,'LipidBlast_negative.xlsx');