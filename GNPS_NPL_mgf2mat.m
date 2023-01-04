clearvars
clc
close all
fclose all;
instrument_type={'LC-ESI-MAXIS';
    'LC-ESI-QTOF';
    'LC-ESI-ORBITRAP'};
datasets={
    'GNPS_NPL_MAXIS';
    'GNPS_NPL_QTOF';
    'GNPS_NPL_ORBITRAP'};
% ---------------------------------------------------------------------
load test spectra
mgfStruct = readNatProdMGF('GNPS-NIH-NATURALPRODUCTSLIBRARY.mgf');
R1PTable=struct2table(mgfStruct.scan);
mgfStruct = readGNPS_NatProdR2_MGF('GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_NEGATIVE.mgf');
R2NTable=struct2table(mgfStruct.scan);
mgfStruct = readGNPS_NatProdR2_MGF('GNPS-NIH-NATURALPRODUCTSLIBRARY_ROUND2_POSITIVE.mgf');
R2PTable=struct2table(mgfStruct.scan);
GNPSTable=[R1PTable;R2NTable;R2PTable]; % table items will be aligned automatically
% save GNPS_NIH_NatureProducts GNPSTable
% load GNPS_NIH_NatureProducts
for i=1:length(instrument_type)
    useid=(strcmpi(GNPSTable.adduct,'M-H') | strcmpi(GNPSTable.adduct,'M+H') | strcmpi(GNPSTable.adduct,'[M+H]+')) ...
        & GNPSTable.isCHONSP & (GNPSTable.monisotopic_molecular_weight <= 1500) & strcmpi(GNPSTable.instrument_type,instrument_type{i});
    GNPSTable_used=GNPSTable(useid,:);
    spectnum=size(GNPSTable_used,1); % total number of spectra
    spectra=table2struct(GNPSTable_used);
    save([datasets{i},'_spectra'],'spectra');
end

