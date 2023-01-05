function mgfStruct = readGNPS_NatProdR2_MGF(file)
% INPUT
% A character array (1 file) or cell array (multiple files) are optional
% input. When the file path is not given as input, files can be browsed for
% 
% OUTPUT
% The output is a structure containing the following fields:
% - scanName: Title given by ProteoWizard (character array)
% - precursorMass: Mass of precursor ion (double)
% - precursorIntensity: Intensity of precursor ion (double)
% - z: Charge of precursor ion (double)
% - scanData: MS/MS spectrum (N x 2 cell containing m/z values
% and intensities)
% 
% Cite As
% Joris Meurs (2022). A MATLAB reader for MASCOT Generic Format (.mgf) files (https://www.mathworks.com/matlabcentral/fileexchange/76033-a-matlab-reader-for-mascot-generic-format-mgf-files), MATLAB Central File Exchange. Retrieved June 13, 2022.
% Modified by
% Author: Ke-Shiuan Lynn Ph.D.
% Assistant Professor
% Department of Mathematics
% Fu-Jen Catholic University
% Email: 128171@mail.fju.edu.tw
% Final Update: Jul. 24, 2022


% Validate file
ismgf(file);
% ------------------ modification start ----------------------------
% allowed elements
KeepElem={'C','H','Br','Cl','F','I','N','O','P','S','Si'};
% ------------------ modification end ----------------------------
fprintf('Reading file\n');
% Open and read .MGF file
fileID = fopen(file,'r');   
fileData = textscan(fileID,'%s');
fileData = fileData{1,1};
fclose(fileID);
% Get total number of scans
scanCount = numel(find(contains(fileData,'BEGIN')));
fprintf('Number of MS/MS scans: %d\n',scanCount);
% Get scan names
fprintf('Obtaining scan titles...\n');
try
    titleRow = find(contains(fileData,'TITLE='));
    for n = 1:length(titleRow)
       tempTitle = fileData{titleRow(n),1};
       titleIDX = find(tempTitle=='=');
       mgfStruct.scan(n).scanName = tempTitle(titleIDX+1:end);
    end
catch
    warning('Titles not found');
end
% Get precursor mass & intensity
fprintf('Obtaining precursor data...\n');
PrecursorRow = find(contains(fileData,'PEPMASS='));
% ------------------ modification start ----------------------------
OrganismRow = find(contains(fileData,'ORGANISM='));
PIRow = find(contains(fileData,'PI='));
IonModeRow = find(contains(fileData,'IONMODE='));
SMILESRow = find(contains(fileData,'SMILES='));
FormulaRow = find(contains(fileData,'INCHI='));
AssessionRow = find(contains(fileData,'SPECTRUMID='));
InstructRow = find(contains(fileData,'SOURCE_INSTRUMENT='));
% ------------------ modification end ----------------------------
try
    for n = 1:length(PrecursorRow)
        tempPrecursor = fileData{PrecursorRow(n),1};
        precursorIDX = find(tempPrecursor=='=');
        mgfStruct.scan(n).precursorMass = double(str2double(tempPrecursor(precursorIDX+1:end)));
        % Get Name string
        tempName = strjoin(fileData((OrganismRow(n)+1):(PIRow(n)-1)),' ');
        tempName = strtrim(erase(tempName,'"'));
        % extract adduct info
        spid=strfind(tempName,' ');
        tempadduct=tempName(spid(end)+1:end); % extract the last sentence sepatated by a space
        if contains(tempadduct,'M+') || contains(tempadduct,'M-')
            if contains(tempadduct,'[') && contains(tempadduct,']') 
                id1=strfind(tempadduct,'[');
                id2=strfind(tempadduct,']');
                mgfStruct.scan(n).adduct = tempadduct((id1+1):(id2-1));
            else
                mgfStruct.scan(n).adduct = tempadduct;
            end
            tempName=tempName(1:spid(end)-1);
        else
            mgfStruct.scan(n).adduct = 'unknown';
        end
        % extract compound name
        if contains(tempName,'!')
            tempstr=strsplit(strtrim(tempName),'!');
            if ~isempty(tempstr{2})
                mgfStruct.scan(n).name = tempstr{2};
            else
                mgfStruct.scan(n).name = tempstr{1}(6:end);
            end
        else
            mgfStruct.scan(n).name = tempName;
        end
        rid=strfind(mgfStruct.scan(n).name,'[IIN-based');
        if ~isempty(rid) && (rid > 5)
            mgfStruct.scan(n).name=strtrim(mgfStruct.scan(n).name(1:(rid-1)));
        end
        tempname=mgfStruct.scan(n).name;
        if contains(tempname,'NAME=')
            id=strfind(tempname,'_');
            if length(id) < 2
                id1=strfind(tempname,' (');
                if isempty(id1)
                    tempname=strrep(tempname,'NAME=','');
                else
                    tempname=tempname(6:id1-1);
                end
            else
                if isempty(tempname((id(end)+1):end))
                    tempname=tempname(6:id(1)-1);
                else
                    tempname=tempname((id(end)+1):end);
                end
            end
            mgfStruct.scan(n).name=tempname;
        end
%         if isempty(mgfStruct.scan(n).name)
%             disp('here')
%         end
        % Get ion mode
        tempIonMode = fileData{IonModeRow(n),1};
        IonModeIDX = find(tempIonMode=='=');
        mgfStruct.scan(n).ionization_mode=lower(tempIonMode(IonModeIDX+1:end));
        if strcmp(mgfStruct.scan(n).adduct,'M+H')
            mgfStruct.scan(n).monisotopic_molecular_weight=mgfStruct.scan(n).precursorMass-1.007276;
        elseif strcmp(mgfStruct.scan(n).adduct,'M-H')
            mgfStruct.scan(n).monisotopic_molecular_weight=mgfStruct.scan(n).precursorMass+1.007276;
        else
            mgfStruct.scan(n).monisotopic_molecular_weight=-1;
        end
        % Get Instrument Type
        tempInstrument = fileData{InstructRow(n),1};
        InstrumentIDX = find(tempInstrument=='=');
        mgfStruct.scan(n).instrument_type=upper(tempInstrument(InstrumentIDX+1:end));
        % Get SMILES and convert to chemical formula
        tempSMILES = fileData{SMILESRow(n),1};
        SMILESIDX = find(tempSMILES=='=');
        SMILES=tempSMILES(SMILESIDX(1)+1:end);
        mgfStruct.scan(n).SMILES=SMILES;
        % Extract chemical formula from INCHI 
        tempFormula = fileData{FormulaRow(n),1};
        tempstr=strsplit(tempFormula,'/');
        mgfStruct.scan(n).chemical_formula=tempstr{2};
        EleList=regexp(mgfStruct.scan(n).chemical_formula,['[','A':'Z','][','a':'z',']?'],'match');
        isCHONSP=~contains(EleList,KeepElem([3,4,5,6,11]));
        mgfStruct.scan(n).isCHONSP=all(isCHONSP);
        tempAssession = fileData{AssessionRow(n),1};
        AsssessionIDX = find(tempAssession=='=');
        mgfStruct.scan(n).assession = tempAssession(AsssessionIDX+1:end);
    end
catch 
    warning('Precursor data not found')
end
% Get scan data and charge (if present)
fprintf('Obtaining MS/MS scans...\n');
beginRow = find(contains(fileData,'BEGIN'));
endRow = find(contains(fileData,'END'));
if length(beginRow) ~= length(endRow)
   error('Unequal rows'); 
end
for n = 1:length(beginRow)
    tempScanData = fileData(beginRow(n):endRow(n),1);
    chargeRow = find(contains(tempScanData,'CHARGE='));
    if ~isempty(chargeRow)
        tempCharge = fileData{chargeRow,1};
        chargeIDX = find(tempCharge=='='); 
        %mgfStruct.scan(n).z = double(str2double(tempCharge(chargeIDX+1:end-1)));
        mgfStruct.scan(n).z = double(str2double(tempCharge(chargeIDX+1:end)));
        parameterRow = find(contains(tempScanData,'='));
        scanRowStart = parameterRow(end)+1;
        tempMS = str2double(tempScanData(scanRowStart:end-1,1));
        tempMZ = tempMS(1:2:end);
        tempINT = tempMS(2:2:end);
        %tempScan = [tempMZ,tempINT];
    else
        mgfStruct.scan(n).z = NaN; 
        scanRowStart = find(contains(tempScanData,'PEPMASS='));
        scanRowStart = scanRowStart+2;
        try
            tempMS = str2double(tempScanData(scanRowStart:end-1,1));
            tempMZ = tempMS(1:2:end);
            tempINT = tempMS(2:2:end);
            %tempScan = [tempMZ,tempINT];
        catch
            tempMS = str2double(tempScanData(scanRowStart-1:end-1,1));
            tempMZ = tempMS(1:2:end);
            tempINT = tempMS(2:2:end);
            %tempScan = [tempMZ,tempINT];
        end
    end
    %mgfStruct.scan(n).scanData = tempScan;
    mgfStruct.scan(n).peak_counter=length(tempMZ);
    mgfStruct.scan(n).ms2_peak=['0;',num2str(tempINT(1)),';',num2str(tempMZ(1)),';0'];
    for i=2:mgfStruct.scan(n).peak_counter
        mgfStruct.scan(n).ms2_peak=strcat(mgfStruct.scan(n).ms2_peak,';',['0;',num2str(tempINT(i)),';',num2str(tempMZ(i)),';0']);
    end
end
fprintf('Finished!\n');
end
function ismgf(files)   
    if ~iscell(files)
        tf = ~isempty(find(contains(files,'mgf'), 1));
        if tf == true
            fprintf('File type: .mgf\n');
        else
            error('File type is not .mgf');
        end
    else
        for p = 1:length(files)
            tf = ~isempty(find(contains(files{p},'mgf'), 1));
            if tf == true
                fprintf('File %d. Type is .mgf\n',p);
            else
                error('File %d. Type is not .mgf',p);
            end
        end
    end  
end
