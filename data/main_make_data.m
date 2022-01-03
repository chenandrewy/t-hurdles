% convert Chen Zim csvs and xlsx to matlab 
% 2021 11 ac, cleaned up previosu stuff

clear

%% === create cz summary data ===
% had to save xlsx as xls to get matlab to read in
% = begin with header info, keep only predictors

% read
temphead = readtable('CZ_data_2021_04/SignalDocumentation.xls','Sheet','BasicInfo');
temphead2 = readtable('CZ_data_2021_04/SignalDocumentation.xls','Sheet','AddInfo');
tempsum = readtable('CZ_data_2021_04/PredictorSummary.xlsx');
tempsum = tempsum(:,{'signalname','tstat'});
tempsum.Properties.VariableNames{1} = 'Acronym';

% merge
tempmerge = innerjoin(temphead, temphead2);
tempmerge = innerjoin(tempmerge, tempsum);

% keep only good reproductions
iok = tempmerge.tstat > 1.5;
tempmerge = tempmerge(iok,:);

signalname = tempmerge.Acronym;
samplebegin = tempmerge.SampleStartYear;
sampleend   = tempmerge.SampleEndYear;
pubyear = tempmerge.Year;
predop  = tempmerge.PredictabilityInOP;

czinfo = table(signalname, samplebegin, sampleend, pubyear, predop);


%% === import and clean cz returns ===

% import
raw = readtable('CZ_data_2021_04/PredictorLSretWide.csv');

% remove poor reproductions
repok = raw;
iok = zeros(1,206);
iok(1) = 1;
for i = 2:206
    if any(strcmp(czinfo.signalname, repok.Properties.VariableNames{i}))        
        iok(i) = 1;
    end
end
repok = repok(:,logical(iok));

% fullsamp sample returns
fullsamp = repok;
for i = 2:size(repok,2)
    if ~isnumeric(fullsamp{:,i})
        fullsamp.(i) = str2double(fullsamp{:,i});
    end
end

%% in sample only
% non-insamp are set to nan
insamp = fullsamp;
y = year(fullsamp.date);
clear clearpred
for i = 2:size(repok,2)
    j = strcmp(czinfo.signalname,insamp.Properties.VariableNames{i});      
    tgood = y >= czinfo(j,:).samplebegin & y <= czinfo(j,:).sampleend;
    insamp{~tgood,i} =  nan;         
end

%% calculate moments and their standard errors

% standard errors are inaccurate here due to pub bias, but they should be
% good enough for weighting in SMM

NbootQ = 200;
seed = 1;

rmat = table2array(insamp(:,2:end)); 

[tpub,rpub,sigpub,corrpub,corrpubmat,nobspub] = rmat2pubdat(rmat);
[momt,momr,momsig,momcorr,Npub] = momfun(tpub,rpub,sigpub,corrpub);

rng(seed,'CombRecursive');

Tboot = size(rmat,1);
Nboot = size(rmat,2);

tic


parfor booti = 1:NbootQ
    
    if mod(booti,20)==0
        disp(['booti = ' int2str(booti)])
    end    

    % draw sample        
    imonth = randi(Tboot, [Tboot 1]);
    isignal = randi(Nboot, [Nboot 1]);    
    
    temprmat = rmat(imonth,:);
    temprmat = temprmat(:,isignal);
    
    % moments
    [temptpub,temprpub,tempsigpub,tempcorrpub] = rmat2pubdat(temprmat);
    
    [boott(booti,:),bootr(booti,:),bootsig(booti,:),bootcorr(booti,:)] ...
        = momfun(temptpub,temprpub,tempsigpub,tempcorrpub);


end
toc

%% pack up

est.t = momt;
est.r = momr;
est.sig = momsig;
est.corr = momcorr;

se.t = std(boott);
se.r = std(bootr);
se.sig = std(bootsig);
se.corr = std(bootcorr);

clear czmom
czmom.est = est;
czmom.se = se;

czsamp = packstruct(tpub,rpub,sigpub,corrpub,corrpubmat,nobspub);


%% === save to disk ===
save('cz-panel.mat','fullsamp','insamp','czmom','czsamp')



%% check

head(fullsamp)

head(insamp)

size(insamp)