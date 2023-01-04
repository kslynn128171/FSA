clearvars
clc
close all

mgfStruct = readGNPS_PNNL_LIPIDS_MGF('PNNL-LIPIDS-NEGATIVE.mgf');
PNNLNTable=struct2table(mgfStruct.scan);
% save GNPS_PNNL_LIPIDS_Neg PNNLNTable
mgfStruct = readGNPS_PNNL_LIPIDS_MGF('PNNL-LIPIDS-POSITIVE.mgf');
PNNLPTable=struct2table(mgfStruct.scan);
% save GNPS_PNNL_LIPIDS_Pos PNNLPTable
% load GNPS_PNNL_LIPIDS_Neg.mat
% load GNPS_PNNL_LIPIDS_Pos.mat
PNNLTable=[PNNLNTable;PNNLPTable];
useid=(strcmpi(PNNLTable.adduct,'[M-H]-') | strcmpi(PNNLTable.adduct,'[M+H]+')) ...
    & contains(PNNLTable.instrument_type,'HCD') & (PNNLTable.monisotopic_molecular_weight <= 1500);
PNNLTable=PNNLTable(useid,:);
save PNNLLIPID_spectra PNNLTable
