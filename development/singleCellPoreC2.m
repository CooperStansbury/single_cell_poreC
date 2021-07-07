% GLOBAL VARIABLES
INDIR='problem_region/';
INFILE='aligns_example_problem.mat';
INPATH=strcat(INDIR, INFILE);

CSV_OUTDIR='problem_region/';
CSV_OUTFILE='aligns_example_problem.csv';
CSV_OUTPATH=strcat(CSV_OUTDIR, CSV_OUTFILE);

INCIDENCE_OUTDIR='problem_region/';
INCIDENCE_OUTFILE='aligns_example_problem_incidence.csv';
INCIDENCE_OUTPATH=strcat(INCIDENCE_OUTDIR, INCIDENCE_OUTFILE);

CHROMOSOME=19;

% Load the dataset
IN_MAT = load(INPATH);
fields = fieldnames(IN_MAT);
aligns=IN_MAT.(fields{1}); 

% write the dataset as a csv
writetable(aligns, CSV_OUTPATH);

% Create fragmnent mid-point markers
aligns.fragment_mean = ceil((aligns.fragment_start+aligns.fragment_end)/2);

% Extract porecMatrix
porecMatrix = table2array(aligns(aligns.chrom==CHROMOSOME, {'read_idx','chrom', 'fragment_mean'}));

% Remove repeated rows
porecMatrix = unique(porecMatrix, 'rows'); 

% Construct incidence matrix
[~, ~, readIDReIdx] = unique(porecMatrix(:, 1)); % Reindex the ReadID
[uniqueContacts, ~, lociReIdx] = unique(porecMatrix(:, 3)); % Reindex the loci
incidenceMatrix = sparse(lociReIdx, readIDReIdx, 1);
incidenceMatrix = unique(incidenceMatrix', 'rows')'; % Remove repeated edges

% Compute order distribution
edgeOrder = full(sum(incidenceMatrix, 1));
figure, histogram(edgeOrder);

% Extract all hyperedges
highOrderContacts = cell(size(incidenceMatrix, 2), 1);
for i = 1:size(incidenceMatrix, 2)
    highOrderContacts{i} = uniqueContacts(find(incidenceMatrix(:, i)>0)); %#ok<FNDSB> % Recover the original loci 
end
highOrderContacts = flipud(highOrderContacts(cellfun('length', highOrderContacts)>=2));

% Create csv file for PAOHvis
CSVIncidenceMatrix = [];
for i = 1:length(highOrderContacts)
    contact = highOrderContacts{i};
    for j = 1:length(contact)
        CSVIncidenceMatrix = [CSVIncidenceMatrix; [i contact(j)]]; %#ok<AGROW>
    end 
end

label = CHROMOSOME*ones(size(CSVIncidenceMatrix, 1), 1);
CSVIncidenceMatrixTable = table(CSVIncidenceMatrix(:, 1), num2str(CSVIncidenceMatrix(:, 2), '%08d'), ...
    cell(size(CSVIncidenceMatrix, 1), 1), cell(size(CSVIncidenceMatrix, 1), 1), label);
writetable(CSVIncidenceMatrixTable, INCIDENCE_OUTPATH)

