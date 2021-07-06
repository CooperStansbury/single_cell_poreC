clear
clc
close all

% input stuff;
expr='02';                                                                 % folder with parquet files (doubles as runID)
%vars={'read_idx','chrom','start','end','strand','mapping_quality','fragment_start','fragment_end','pass_filter'};

% data stuff
rt=readtable('GRCm39.txt');                                                % reference genome
kys=rt.RefSeq;
vls=rt.Chrom;
n2c=containers.Map(kys,vls);

pds=parquetDatastore(expr);
%pds.SelectedVariableNames=vars;

fnums=pds.Files;
fnums=regexprep(fnums,'^.*_batch','');
fnums=regexprep(fnums,'_.*','');
fnums=double(string(fnums));
[~,idx]=sort(fnums);

sorted_files=pds.Files(idx);

%pds.Files=pds.Files(idx);
%pds=partition(pds,'Files',pds.Files);

runID=double(string(expr));

i=0;
FT=table;
offset=0;
max_idx=50000;
%while(pds.hasdata)
for i=1:length(sorted_files)
    parquet_file=sorted_files{i};
    fprintf("%s\n",parquet_file);
    %T=pds.read;
    T=parquetread(parquet_file);
    %T=parquetread(parquet_file,'SelectedVariableNames',vars);
    abs_max=double(max(T.read_idx));                                       % sigh
    
    T=T(T.pass_filter,1:end-1);
    idx=ismember(T.chrom,kys);
    T=T(idx,:);
    T.chrom=cell2mat(values(n2c,cellstr(T.chrom)));
    T=addvars(T,runID*ones(height(T),1),'Before',1,'NewVariableNames',{'runID'});
    T.Properties.VariableNames{'start'}='xStart';
    
    %{
    for i=1:width(T)
        T.(T.Properties.VariableNames{i})=double(T.(T.Properties.VariableNames{i}));
    end
    %}
    
    % they RESET the read_idx every loop through (WHY!)
  
    %max_idx=max(T.read_idx)+1;
    T.read_idx=T.read_idx+offset;
    offset=offset+max_idx;
    
    FT=[FT;T];                                                             % rename and save
end
disp(i)
