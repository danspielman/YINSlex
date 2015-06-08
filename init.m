% script init.m
%
% run from this directory to set up the paths correctly.

addpath(genpath(fullfile(pwd,'/matlab')));

% initialize our java classes

javaaddpath out/YINSlex_jar/YINSlex.jar

import lexvolts.*;
import lexvoltsdirected.*;