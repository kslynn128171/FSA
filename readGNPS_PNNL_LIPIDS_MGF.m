function mgfStruct = readGNPS_PNNL_LIPIDS_MGF(file)
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
AssessionRow = find(contains(fileData,'SPECTRUMID='));
InstructRow = find(contains(fileData,'SOURCE_INSTRUMENT='));
% ------------------ modification end ----------------------------
try
    for n = 1:length(PrecursorRow)
        % Get ion mode
        tempIonMode = fileData{IonModeRow(n),1};
        IonModeIDX = find(tempIonMode=='=');
        mgfStruct.scan(n).ionization_mode=lower(tempIonMode(IonModeIDX+1:end));
        % Get Name string
        tempName = strjoin(fileData((OrganismRow(n)+1):(PIRow(n)-1)),' ');
        tempName = strtrim(erase(tempName,';'));
        tempstr=strsplit(tempName,' ');
        id=2;
        while id < length(tempstr)
            if ~strcmpi(tempstr{id},'[M+H]+') && ~strcmpi(tempstr{id},'[M-H]-')
                id=id+1;
            else
                break;
            end
        end
        if strcmpi(tempstr{id},'[M+H]+') || strcmpi(tempstr{id},'[M-H]-')
            % extract compound name
            if contains(tempstr{id-1},'NAME')
                mgfStruct.scan(n).name = tempstr{id-1}(6:end);
            else
                mgfStruct.scan(n).name = tempstr{id-1};
            end
            % extract adduct info
            mgfStruct.scan(n).adduct=tempstr{id};
            % Extract chemical formula from SMILES 
            mgfStruct.scan(n).chemical_formula=tempstr{id+1};
            [NumID,~]=regexp(mgfStruct.scan(n).chemical_formula,'H\d*[A-Za-z]','match');
            HNum=str2double(NumID{1}(2:end-1));
            % monisotopic_molecular_weight
            tempPrecursor = fileData{PrecursorRow(n),1};
            precursorIDX = find(tempPrecursor=='=');
            mgfStruct.scan(n).monisotopic_molecular_weight=double(str2double(tempPrecursor(precursorIDX+1:end)));
            % precursorMZ
            if contains(mgfStruct.scan(n).adduct,'M+H') % positive mode
                mgfStruct.scan(n).precursorMass = mgfStruct.scan(n).monisotopic_molecular_weight+1.007276;
                % adjust chemical formula
                tempform=strrep(mgfStruct.scan(n).chemical_formula,NumID{1},['H',num2str(HNum-1),NumID{1}(end)]);
                kid=true(size(tempform));
                [~,NumID]=regexp(tempform,'[A-Za-z]1[A-Za-z]','match');
                kid(NumID+1)=false;
                if isletter(tempform(end-1)) && (tempform(end)=='1')
                    kid(end)=false;
                end
                mgfStruct.scan(n).chemical_formula=tempform(kid);
            elseif contains(mgfStruct.scan(n).adduct,'M-H') % negative mode
                mgfStruct.scan(n).precursorMass = mgfStruct.scan(n).monisotopic_molecular_weight-1.007276;
                % adjust chemical formula
                tempform=strrep(mgfStruct.scan(n).chemical_formula,NumID{1},['H',num2str(HNum+1),NumID{1}(end)]);
                kid=true(size(tempform));
                [~,NumID]=regexp(tempform,'[A-Za-z]1[A-Za-z]','match');
                kid(NumID+1)=false;
                if isletter(tempform(end-1)) && (tempform(end)=='1')
                    kid(end)=false;
                end
                mgfStruct.scan(n).chemical_formula=tempform(kid);
            else
                mgfStruct.scan(n).precursorMass=-1;
            end
        else
            % extract compound name
            if contains(tempstr{1},'NAME')
                mgfStruct.scan(n).name = tempstr{1}(6:end);
            else
                mgfStruct.scan(n).name = tempstr{1};
            end
            % extract adduct info
            mgfStruct.scan(n).adduct=tempstr{2};
            % Extract chemical formula from SMILES 
            mgfStruct.scan(n).chemical_formula=tempstr{3};
            % monisotopic_molecular_weight
            mgfStruct.scan(n).monisotopic_molecular_weight=-1;
            mgfStruct.scan(n).precursorMass=-1;
        end
        % Get Instrument Type
        tempInstrument = fileData{InstructRow(n),1};
        sidx = find(tempInstrument=='=');
        eidx = find(tempInstrument==';');
        mgfStruct.scan(n).instrument_type=upper(tempInstrument((sidx+1):(eidx-1)));
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
