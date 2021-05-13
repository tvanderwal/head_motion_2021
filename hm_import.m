%% hm_import.m
% Simon Frew | NNL | BCCHRI
% create data struct out of MCFLIRT par files and HBN behavioural data. 


% final output: 
%   hm_data with fields: 
%       .id             : sub-XXXXXXXXXX
%       .sex            : From import_hbnpheno()
%                           0=male
%                           1=female
%       .age            : From import_hbnpheno()
%       .site           : Study site. From import_hbnpheno()
%       .CBCL_Total_T   : From import_CBCL() (if available)
%       .rest1          : Raw displacement data from Rest1. Imported from MCFLIRT .par files
%       .rest2          : Raw displacement data from Rest2. Imported from MCFLIRT .par files
%       .movieDM        : Raw displacement data from MovieDM. Imported from MCFLIRT .par files

fprintf("hm_import.m\n")

%% import demographic data
fprintf("\t 1/5\timporting demographic data...\n")
basicdemos = import_basicdemos(); 
CBCL = import_CBCL(); 
hbnpheno = import_hbnpheno(); 


%% import .par motion data to hm_data struct
par_path = fullfile(".", "data", "mcflirt_output"); % path to par files

file_list = dir(par_path);
file_list(1:2) = [];
file_list = file_list([file_list.isdir]==0);

hm_data = struct('id', {}, 'sex', {}, 'age', {}, 'site', {}, 'CBCL_Total_T', {}, 'rest1', {}, 'rest2', {}, 'movieDM', {});

fprintf("\t 2/5\timporting motion data...\n")

for k = 1 : length(file_list) 
    filename = file_list(k).name;
    split_name = strsplit(file_list(k).name, "_");
    % for each file 
    
    % check if subject exists
    % if subject doesnt exist, create
    sub_exists = ismember(string([hm_data.id]), split_name(1));
    

    
    if any(sub_exists)
        
        % check name of file and add par values 
        sub_index = find(sub_exists);
        if contains(split_name(2), "movieDM")
            hm_data(sub_index).movieDM = import_par(fullfile(par_path, filename));
        else
            if contains(split_name(3), "1")
                hm_data(sub_index).rest1 = import_par(fullfile(par_path, filename));
            elseif contains(split_name(3), "2")
                hm_data(sub_index).rest2 = import_par(fullfile(par_path, filename));
            end
        end
        
    else
        hm_data(end+1).id = string(split_name(1));

        % check name of file and add par values 
        if contains(split_name(2), "movieDM")
            hm_data(end).movieDM = import_par(fullfile(par_path, filename));
        else
            if contains(split_name(3), "1")
                hm_data(end).rest1 = import_par(fullfile(par_path, filename));
            elseif contains(split_name(3), "2")
                hm_data(end).rest2 = import_par(fullfile(par_path, filename));
            end
        end
    
    end
end

%% append phenotypic and behavioural data
fprintf("\t 3/5\tadding phenotypic data...\n")

missing_hbnpheno = 0;
for m = 1: length(hm_data)
    sub_id = extractAfter(hm_data(m).id, 4);
    sub_exists = find(ismember(hbnpheno.EID, sub_id), 1, 'first');
    if isempty(sub_exists)
        missing_hbnpheno = missing_hbnpheno + 1;
        continue
    end
    % 0 = male, 1 = female
    hm_data(m).sex = table2array(hbnpheno(sub_exists, 2));
    hm_data(m).age = table2array(hbnpheno(sub_exists, 3));

end 

% walk through hbn and pull CBCL Total T
missing_CBCL = 0;
for m = 1: length(hm_data)
    sub_id = extractAfter(hm_data(m).id, 4);
    sub_exists = ismember(CBCL.EID, sub_id); 
    if ~any(sub_exists)
        missing_CBCL = missing_CBCL + 1;
        continue
    end

    % add to exclude NAN CBCL scores
    cbclVal = table2array(CBCL(sub_exists, end)); 
    if ~isnan(cbclVal)
        hm_data(m).CBCL_Total_T = cbclVal;
    else
        missing_CBCL = missing_CBCL + 1;
        continue
    end

end 

% walk through hbn and pull site
missing_site = 0;
for m = 1: length(hm_data)
    sub_id = extractAfter(hm_data(m).id, 4);
    sub_exists = ismember(basicdemos.EID, sub_id); 
    if ~any(sub_exists)
        missing_site = missing_site + 1;
        continue
    end
    
    hm_data(m).site = basicdemos.StudySite(sub_exists); 

end 

%% remove subjects with missing data 
fprintf("\t 4/5\tremoving subjects with missing values...\n")

missing_sex = 0;
missing_age = 0;
missing_CBCL = 0;
for m = 1: length(hm_data)
    missing_sex = missing_sex + isempty(hm_data(m).sex);
    missing_age = missing_age + isempty(hm_data(m).age);
    missing_CBCL = missing_CBCL + isempty(hm_data(m).CBCL_Total_T);
end

% remove subjects with missing age scores
emptyIndex = find(arrayfun(@(hm_data) isempty(hm_data.age), hm_data)); 
hm_data(emptyIndex) = []; % should be 1 x 1977 struct

% remove subjects missing any raw data...
rest1Empty = find(arrayfun(@(hm_data) isempty(hm_data.rest1), hm_data));
hm_data(rest1Empty) = [];
rest2Empty = find(arrayfun(@(hm_data) isempty(hm_data.rest2), hm_data));
hm_data(rest2Empty) = [];
movieDMEmpty = find(arrayfun(@(hm_data) isempty(hm_data.movieDM), hm_data));
hm_data(movieDMEmpty) = [];

% remove subjects with missing or excess run data
rest1wV = find(arrayfun(@(hm_data) any(size(hm_data.rest1) ~= [375, 6]), hm_data));
hm_data(rest1wV) = [];
rest2wV = find(arrayfun(@(hm_data) any(size(hm_data.rest2) ~= [375, 6]), hm_data));
hm_data(rest2wV) = [];
movieDMwV = find(arrayfun(@(hm_data) any(size(hm_data.movieDM) ~= [750, 6]), hm_data));
hm_data(movieDMwV) = [];

rest1_wV = 0;
rest2_wV = 0;
movieDM_wV = 0;

for i = 1:length(hm_data)
    if any(size(hm_data(i).rest1) ~= [375, 6])
        rest1_wV = rest1_wV + 1; 
    end   
    if any(size(hm_data(i).rest2) ~= [375, 6])
        rest2_wV = rest2_wV + 1; 
    end
    if any(size(hm_data(i).movieDM) ~= [750, 6])
        movieDM_wV = movieDM_wV + 1;
    end

end


%% cleanup
fprintf("\t 5/5\tsaving hm_data.mat\n")

clearvars -except hm_data

save("hm_data.mat")
%% helper functions: 

function par_table = import_par(filename, dataLines)
    %IMPORTFILE Import data from a text file
    %  par_table = import_par(FILENAME) reads data
    %  from text file FILENAME for the default selection.  Returns the data
    %  as a table.
    %
    %  par_table = import_par(FILE, DATALINES) reads
    %  data for the specified row interval(s) of text file FILENAME. Specify
    %  DATALINES as a positive scalar integer or a N-by-2 array of positive
    %  scalar integers for dis-contiguous row intervals.
    %
    %  Example:
    %  par_table = import_par("~\hm_data\par_files\sub-XXXXXXXX_task-XXXXX_bold_mcf.par", [1, Inf]);
    %
    %  See also READTABLE.
    %
    % Auto-generated by MATLAB on 06-Jan-2020 10:25:37

    %% Input handling

    % If dataLines is not specified, define defaults
    if nargin < 2
        dataLines = [1, Inf];
    end

    %% Setup the Import Options
    opts = delimitedTextImportOptions("NumVariables", 6);

    % Specify range and delimiter
    opts.DataLines = dataLines;
    opts.Delimiter = "  ";

    % Specify column names and types
    opts.VariableNames = ["x_rad", "y_rad", "z_rad", "x_mm", "y_mm", "z_mm"];
    opts.VariableTypes = ["double", "double", "double", "double", "double", "double"];
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    par_table = readtable(filename, opts);

end

function basicdemos = import_basicdemos()
    %% Import data from text file
    % Note: This script has been modified in this relesase to preserve privacy as per our data usage aggreement
    % Script for importing data from the following text file:
    %
    %    filename: fullfile(".", "data", "demographics", "XXXXXXXX.csv")v
    %
    % Auto-generated by MATLAB on 06-Jan-2020 12:46:55. Modified 2021/04/19

    filename = fullfile(".", "data", "demographics", "XXXXXXXX.csv");
    %% Setup the Import Options
    opts = delimitedTextImportOptions("NumVariables", 15);

    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["AnonymizedID", "SubjectType", "Visit", "Dayssinceenrollment", "EID", "START_DATE", "Patient_ID", "EnrollmentYear", "EnrollmentSeason", "Sex", "Age", "StudySite", "CommercialUse", "ReleaseNumber", "Participant_Status"];
    opts.VariableTypes = ["double", "categorical", "double", "double", "string", "categorical", "string", "double", "categorical", "double", "double", "double", "categorical", "double", "categorical"];
    opts = setvaropts(opts, [5, 7], "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 3], "TrimNonNumeric", true);
    opts = setvaropts(opts, [1, 3], "ThousandsSeparator", ",");
    opts = setvaropts(opts, [2, 5, 6, 7, 9, 13, 15], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    basicdemos = readtable(filename, opts);

end

function CBCL = import_CBCL()
    %% Import data from text file
    % Note: This script has been modified in this relesase to preserve privacy as per our data usage aggreement 
    % Script for importing data from the following text file:
    %
    %    filename: fullfile(".", "data", "demographics", "XXXXXXXX.csv");
    %
    % Auto-generated by MATLAB on 06-Jan-2020 12:48:48. Modified 2021/04/19

    filename = fullfile(".", "data", "demographics", "XXXXXXXX.csv");

    %% Setup the Import Options
    opts = delimitedTextImportOptions("NumVariables", 157);

    % Specify range and delimiter
    opts.DataLines = [3, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["AnonymizedID", "SubjectType", "Visit", "Dayssinceenrollment", "EID", "Start_Date", "Study", "Site", "Days_Baseline", "Year", "Season", "CBCL_01", "CBCL_02", "CBCL_03", "CBCL_04", "CBCL_05", "CBCL_06", "CBCL_07", "CBCL_08", "CBCL_09", "CBCL_10", "CBCL_11", "CBCL_12", "CBCL_13", "CBCL_14", "CBCL_15", "CBCL_16", "CBCL_17", "CBCL_18", "CBCL_19", "CBCL_20", "CBCL_21", "CBCL_22", "CBCL_23", "CBCL_24", "CBCL_25", "CBCL_26", "CBCL_27", "CBCL_28", "CBCL_29", "CBCL_30", "CBCL_31", "CBCL_32", "CBCL_33", "CBCL_34", "CBCL_35", "CBCL_36", "CBCL_37", "CBCL_38", "CBCL_39", "CBCL_40", "CBCL_41", "CBCL_42", "CBCL_43", "CBCL_44", "CBCL_45", "CBCL_46", "CBCL_47", "CBCL_48", "CBCL_49", "CBCL_50", "CBCL_51", "CBCL_52", "CBCL_53", "CBCL_54", "CBCL_55", "CBCL_56A", "CBCL_56B", "CBCL_56C", "CBCL_56D", "CBCL_56E", "CBCL_56F", "CBCL_56G", "CBCL_56H", "CBCL_57", "CBCL_58", "CBCL_59", "CBCL_60", "CBCL_61", "CBCL_62", "CBCL_63", "CBCL_64", "CBCL_65", "CBCL_66", "CBCL_67", "CBCL_68", "CBCL_69", "CBCL_70", "CBCL_71", "CBCL_72", "CBCL_73", "CBCL_74", "CBCL_75", "CBCL_76", "CBCL_77", "CBCL_78", "CBCL_79", "CBCL_80", "CBCL_81", "CBCL_82", "CBCL_83", "CBCL_84", "CBCL_85", "CBCL_86", "CBCL_87", "CBCL_88", "CBCL_89", "CBCL_90", "CBCL_91", "CBCL_92", "CBCL_93", "CBCL_94", "CBCL_95", "CBCL_96", "CBCL_97", "CBCL_98", "CBCL_99", "CBCL_100", "CBCL_101", "CBCL_102", "CBCL_103", "CBCL_104", "CBCL_105", "CBCL_106", "CBCL_107", "CBCL_108", "CBCL_109", "CBCL_110", "CBCL_111", "CBCL_112", "CBCL_113A", "CBCL_113B", "CBCL_113C", "CBCL_AD", "CBCL_AD_T", "CBCL_WD", "CBCL_WD_T", "CBCL_SC", "CBCL_SC_T", "CBCL_SP", "CBCL_SP_T", "CBCL_TP", "CBCL_TP_T", "CBCL_AP", "CBCL_AP_T", "CBCL_RBB", "CBCL_RBB_T", "CBCL_AB", "CBCL_AB_T", "CBCL_OP", "CBCL_Int", "CBCL_Int_T", "CBCL_Ext", "CBCL_Ext_T", "CBCL_C", "CBCL_Total", "CBCL_TOTAL_T"];
    opts.VariableTypes = ["double", "categorical", "double", "double", "string", "categorical", "categorical", "double", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
    opts = setvaropts(opts, 5, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 3], "TrimNonNumeric", true);
    opts = setvaropts(opts, [1, 3], "ThousandsSeparator", ",");
    opts = setvaropts(opts, [2, 5, 6, 7, 11], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    CBCL = readtable(filename, opts);

end

function hbnpheno = import_hbnpheno()
    %% Import data from text file
    % Note: This script has been modified in this relesase to preserve privacy as per our data usage aggreement    
    % Script for importing data from the following text file:
    %
    %    filename = fullfile(".", "data", "demographics", "XXXXXXXX.csv");
    %
    % Auto-generated by MATLAB on 06-Jan-2020 14:09:52. Modified 2021/04/19

    filename = fullfile(".", "data", "demographics", "XXXXXXXX.csv");

    %% Setup the Import Options
    opts = delimitedTextImportOptions("NumVariables", 6);

    % Specify range and delimiter
    opts.DataLines = [2, Inf];
    opts.Delimiter = ",";

    % Specify column names and types
    opts.VariableNames = ["EID", "Sex", "Age", "EHQ_Total", "Commercial_Use", "Full_Pheno"];
    opts.VariableTypes = ["string", "double", "double", "double", "categorical", "categorical"];
    opts = setvaropts(opts, 1, "WhitespaceRule", "preserve");
    opts = setvaropts(opts, [1, 5, 6], "EmptyFieldRule", "auto");
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";

    % Import the data
    hbnpheno = readtable(filename, opts);
end