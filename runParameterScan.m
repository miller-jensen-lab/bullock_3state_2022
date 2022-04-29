%% RUNPARAMETERSCAN
%{
  This script, when run, will attempt to perform a parameter scan on an
  NFsim model for the selected parameters.  To run this script, update
  the scan options and setup below accordingly to point to the desired
  model and scan the given parameters.  Output of this file is stored in
  the trajectories that are saved to the given output directory.  For an
  example for loop that reads all the output trajectories, see the very
  end of this script.

 Original Code from NFSIM 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   %
%     @@    @  @@@@@                %
%     @ @   @  @                    %
%     @  @  @  @@@@  ___            %
%     @   @ @  @    /__  | |\ /|    %
%     @    @@  @    ___\ | | v |    %
%                                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NFsim - the network free stochastic simulator, v1.11

michael w. sneddon
justin s. hogg
james r. faeder
thierry emonet

Yale University
University of Pittsburgh
funded by the National Science Foundation


################################################################################

NFsim is a free, open-source, biochemical reaction simulator designed to handle systems 
that have a large or even infinite number of possible molecular interactions or states. 
NFsim also has advanced and flexible options for simulating coarse-grained representations 
of complex nonlinear reaction mechanisms.

NFsim is ideal for modeling polymerization, aggregation, and cooperative reactions that 
cannot be handled with traditional stochastic or ODE simulators. Models are specified in 
the BioNetGen Langauge, providing a powerful model building environment.

If you just want to download and use NFsim, you should simply download a preconfigured
packaged release from http://emonet.biology.yale.edu/nfsim.  If you want to hack on the
code or make contributions, please create a fork and submit pull requests to the dev
branch.

If you use NFsim for your research or work, please cite NFsim as:
Sneddon MW, Faeder JR & Emonet T.  Efficient modeling, simulation and 
coarse-graining of biological complexity with NFsim.  Nature Methods,(2011) 8(2):177-83.

################################################################################

Further developed by Margaret Elise Bullock - 04/29/2022
 Running instructions: 
    Edit the .bngl file
    Setup pathway and parameters of interest
    Run this file (which will call runNFsimOnce.m)
    A .mat file will be generated with variables of interest over time
    This can then be analyzed/visualized with the visualization.m file

%}
%%

%  BASIC SETUP.  Change these options below so that they correctly point to
%  the correct directories, models, NFsim installation, and desired output
%  directory.

pathToModel =  'TranscriptionalCyclingModel/';
bnglFileName = 'TranscriptionalCyclingModel.bngl';
pathToNFsim = %'C:\Users\funny\Documents\RuleBender-2.3.2\RuleBender\BioNetGen';
pathToOutput = [pathToModel,'output/'];

%%
%  SCAN OPTIONS.  Add the name of the parameters you would like to scan,
%  togehter with the range and interval of values you would like to scan
%  below.

parameterx = [1 5 10 50 100]; %PBR
parametery = [0.005 0.01 0.05 0.1 .5 1 5]; %BIR
parameterz = [0.005 0.01 0.05 0.1 .5 1 5 10]; %BTR
parameterm = [1 5 10 50 100]; %PPRR

paramNames = {'PBR',  'BIR','BTR','PPRR'};

paramValues = { parameterx; parametery; parameterz; parameterm}; 
            
NumOfRuns=900;

% Information from the .bngl file. Helpful for analysis
endtime = 24; %in hours
stepsize= 1;
savedItems={'Time','UP','AP','BP','RNA','Protein'}; %Same order as Observables in .bngl file

% Save this parameters for access later
save('important','paramNames','stepsize','paramValues','parameterx',...
    'parametery','parameterz','parameterm','NumOfRuns','savedItems','endtime')


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO NOT EDIT THIS SCRIPT BELOW THIS LINE UNLESS YOU UNDERSTAND MATLAB PROGRAMMING

tic;
fprintf('\n+ + + + + + + + + + + + + + + + + + + + + + + + + +\n');
fprintf(['Running parameter scan on file: ',bnglFileName,'\n']);
totalRunCount = 1; currentRunNumber = 1;
for k=1:length(paramValues), totalRunCount = length(paramValues{k})*totalRunCount; end;
fprintf(['Total number of runs that will be generated: ', num2str(totalRunCount),'\n\n' ]);

currentPosition = 1;
currentIndex = ones(size(paramNames));

% run the first trial
fprintf(['\n------\n+ run (',num2str(currentRunNumber),' of ',num2str(totalRunCount),')\n']);
fprintf(['+ running with parameters:\n     ']);
for k=1:length(paramNames),
    fprintf(['   ',paramNames{k}]);
    fprintf(['=',num2str(paramValues{k}(currentIndex(k)))]);
end;
fprintf('\n');
fprintf(['+ trajectories will be saved to file with suffix: ',num2str(currentRunNumber),'.gdat\n']);
paramArray = size(paramNames);
for k=1:length(paramNames), paramArray(k) = paramValues{k}(currentIndex(k)); end;
console = runNFsimOnce(pathToModel,bnglFileName, pathToNFsim, pathToOutput, ...
    currentRunNumber, paramNames, paramArray,NumOfRuns,savedItems,stepsize*endtime+1);
thisRunTime = toc;
meanRunTime = thisRunTime/currentRunNumber;
fprintf(['+ done.  this run took ',num2str(thisRunTime),'s to complete.\n']);
fprintf(['+ estimated time to completion: ',num2str(meanRunTime.*(totalRunCount-currentRunNumber)./60),' min\n']);



while true,
    
     %reset the value at a particular position if we count too high
    allDone = 0;
    while currentIndex(currentPosition) == length(paramValues{currentPosition})
        currentIndex(currentPosition) = 1;
        currentPosition = currentPosition+1;
        %check if we are completely done (once the last digit roles over)
        if currentPosition>length(paramNames), allDone=1; break;  end;
        continue;
    end;
    if allDone, break; end;
    
     %if we got here, then we can increment the value at the current position
    currentIndex(currentPosition) = currentIndex(currentPosition)+1;
    currentPosition = 1;
   
    
    % run the actual simulation
    fprintf(['\n------\n+ run (',num2str(currentRunNumber+1),' of ',num2str(totalRunCount),')\n']);
    fprintf(['+ running with parameters:\n     ']);
    for k=1:length(paramNames), 
        fprintf(['   ',paramNames{k}]);
        fprintf(['=',num2str(paramValues{k}(currentIndex(k)))]);
    end;
    fprintf('\n');
    fprintf(['+ trajectories will be saved to file with suffix: ',num2str(currentRunNumber+1),'.gdat\n']);
    paramArray = size(paramNames);
    for k=1:length(paramNames), paramArray(k) = paramValues{k}(currentIndex(k)); end;
    console = runNFsimOnce(pathToModel,bnglFileName, pathToNFsim, pathToOutput, ...
        currentRunNumber+1, paramNames, paramArray, NumOfRuns,savedItems,stepsize*endtime+1);
    thisRunTime = toc;
    meanRunTime = thisRunTime/currentRunNumber;
    fprintf(['+ done.  this run took ',num2str(thisRunTime),'s to complete.\n']);
    fprintf(['+ estimated time to completion: ',num2str(meanRunTime.*(totalRunCount-currentRunNumber)./60),' min\n']);
    
    
    
    currentRunNumber = currentRunNumber +1;
    k=k+1;
end;
    
    

totalTime = toc;
fprintf(['\n\nParameter Scan Complete.  \nTotal Run time: ',num2str(totalTime./60),' minutes\n']);
fprintf(['Average Time Per Run: ',num2str(meanRunTime./60),' minutes\n']);
fprintf('+ + + + + + + + + + + + + + + + + + + + + + + + + +\n');
