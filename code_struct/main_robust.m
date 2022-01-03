% Andrew 2019 04


restoredefaultpath

% add my tools to path
mytools = genpath(['Matlab_Tools/']);
path(mytools,path);
clear mytools;

clear


%% import results

robustlist = {
    'rhoconst_300_48'
    'rho20_300_48'
    'rhodisp_500_96'
    'riskprem_alt_500_96_300'    
    'logit_500_96'
    'gamma2_500_96'
    'gamma05_500_95_500'
    't2group_500_96_300'
    't2groupalt_500_96_300'
    'sigfat_500_96_500'
};



% import all
disp('importing lots of bootstraps')
for i = 1:length(robustlist)
    name = robustlist{i};
    temp = load(['../output/boot_' name '/results/check_boot_' name '.mat']);
    boot{i} = temp.boot;
end

%% summarize

vlist = {'t_fdr05','t_fdr01','pubfdr','smooth_s'};

% get numbers
clear stat
for i = 1:length(boot)
    tempboot = boot{i};
    for vi = 1:length(vlist)
        v = vlist{vi};
        tempvar = tempboot.(v);
        
        % remove in rare cases when smooth_s is negative
        if strcmp(v,'smooth_s')
            tempvar(tempvar<0) = [];
        end
        
        stat.(v)(i,:) = [
            median(tempvar) std(tempvar)
        ];
    end
    
end

% organize into table
num  = [];
top = {};
for vi = 1:length(vlist)
    num  = [
        num  stat.(vlist{vi})
        ];
    top = [
        top vlist{vi} {''}
        ];
end

tab = [
    top
    num2cell(num )
    ];

tab = [
    [{''}; robustlist] tab
    ]

open tab