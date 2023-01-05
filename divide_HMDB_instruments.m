clearvars
clc
close all
instrument={'QTOF','ITFT','QFT'};
HMDBTable=readtable('HMDB_CHONPS_MSMS.txt','FileType','delimitedtext',...
   'Delimiter','\t','TextType','string','VariableNamingRule','modify');
for i=1:length(instrument)
    Table=HMDBTable(contains(upper(HMDBTable.instrument_type),instrument{i}),:);
    spectra=table2struct(Table);
    save(['HMDB_',instrument{i},'_spectra.mat'],'spectra');
end
