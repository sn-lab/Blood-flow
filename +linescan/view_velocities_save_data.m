         
% view_velocities.m
% goes through files generated by pic_to_vel_dia (either velocity or
% diameter is ok, although all graphs are labeled for velocity)
% Put all files in a folder and the program will loop through all the files
% in that folder.
% add a histogram
% Plot velocity data from pic_to_veldatanewsope and allows user to select
% the range of good velocities
% 12-09-04 Adding Flux capabilities
% 12-23-04 Add stim-triggered abilites. only 1 stim channel, 1 analog
% channel
% Check and compares with official version of view-velocities
% 02-19-05 closing files to deal with memory problems
% 07-26-05 Stim Trig avg
% 09-16-05: Make official version
% 10-20-05: stimthresh for stimulus has to be smaller for some files
% 10-21-05: change the method for pausing, and going throughmany movies
% 10-01-07: changed file name
% 03-11-08: also calculate abs(vels)
% 08-10-11: saves "good" timepoints and velocities in a new file with
% variable names time2 and vel2

display('version 08-10-11')

clear;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Parameters

runavgtime = 5000;
select_range = 1;
Twin = .5; % length of time to look for min/max (s)
threshforstim = 0; % or 50
threshforstim = 100; % FOR 08-10-11 Schaffer 

% Running parameters
prompt={'Go through multiple files in folder?', 'Display fluxes?', 'Stimulus (none, a, b)?', 'Analog input (none, a, b)?'};
def={'yes', 'no', 'none', 'none'};
dlgTitle='Processing parameters';
lineNo=1;
answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');

if strcmp(answer(1), 'yes'); use_many = 1;
else; use_many = 0; end;

if strcmp(answer(2), 'yes'); use_flux = 1;
else; use_flux = 0; end;

if strcmp(answer(3), 'none'); usestim = 0;
elseif strcmp(answer(3),'a') | strcmp(answer(3), 'A');
    usestim=7;
elseif strcmp(answer(3),'b') | strcmp(answer(3), 'B');
    usestim=8;
else usestim = 0;
end

if strcmp(answer(4), 'none'); useana = 0;
elseif strcmp(answer(4),'a') | strcmp(answer(3), 'A');
    useana=7;
elseif strcmp(answer(4),'b') | strcmp(answer(3), 'B');
    useana=8;
else useana = 0;
end

if use_many % write data
    [filename, pathname] = uiputfile('*.csv', 'Comma delimited file save As');
    Datafile = [pathname, filename];
    Datafiletemp = [pathname, 'temp', filename];    
end;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if usestim ~= 0 % get parameters for stim-triggered avg
    prompt={'Tcycle (s)',... 
            'Tbin (s)',... 
            'Tstep (s)'}
    def={'20', '0.5', '0.1'};
    dlgTitle='Plot parameters';
    lineNo=1;
    answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');
    Tcycle = str2double(answer(1));  % 2* minimum cycle length
    Tbin = str2double(answer(2));      % bin to average
    Tstep = str2double(answer(3));
    boxsize = 10;
end;


if use_many
    [fname,pname] = uigetfile('*.*', 'Select a file in the folder of velocity files');
    files = dir([pname '*']);
    
    clearname = [];
    for temp = 1:length(files)
        if files(temp).name(1) == '.'
            clearname = [clearname temp]
        end
    end
    files(clearname) = [];
        
    nmovies = length(files);
    cd(pname);
    
    Names = 'Names';
    Velocities = 'Velocity';
    StdDev = 'Std. Dev.';
    AbsVelocities = 'Abs Vel';
    AbsStdDev = 'Abs Std. Dev.';
    MinRanges = 'Minimum';
    MaxRanges = 'Maximum';
    PerValids = '% Points Valid';
    AvgSyss = 'Systolic Vel.';
    AvgDias = 'Diastolic Vel.';
    if use_flux == 1
        AvgFlux1s = 'Flux1 avg';
        AvgFlux2s = 'Flux2 avg'
        AvgFlux3s = 'Flux3 avg';
    end
    
else
    [fname,pname] = uigetfile('*.*');
    Openfile = [pname, fname];
    load(Openfile);
    nmovies = 1;
end



for i = 1:nmovies % loop through movies
    if use_many
        Openfile = [pname files(i).name];
    end
    disp(['Loading ' Openfile]);
    load(Openfile);
    
    if not(exist('Tfactor'))
        Tfactor = 1.333; % for verH
        %                 Tfactor = 1.43; %verC
    end
    if not(exist('WinPixelsDown'))
        WinPixDown = 32;
    else
        WinPixDown = WinPixelsDown;
        clear WinPixelsDown;
    end
    
    % assign data variables
    time = squeeze(Result(:,2))/1000;
    lines = squeeze(Result(:,1));
    vel = squeeze(Result(:,3));
    if use_flux == 1
        flux1 = squeeze(Result(:,6));
        flux2 = squeeze(Result(:,9));
        flux3 = squeeze(Result(:,10));
    end; % if use_flux
    if usestim~=0 
        timestim = squeeze(Result(:,2))/1000;
        stim = squeeze(Result(:,usestim));
    end; % usestim
    
    % smooth velcocities
    numtoavg = ceil(runavgtime/(WinPixDown/Tfactor));
    convbox = ones([1 numtoavg])'/numtoavg;      
    smoothed = convn(vel, convbox, 'same');
    
    avg = mean(vel);
    stdev = std(vel);
    if use_flux
        avgflux1 = mean(flux1);avgflux2 = mean(flux2);avgflux3 = mean(flux3);
    end
    
    Y = fft(vel, 512);
    Pyy = Y.* conj(Y) / 512;
    f = 1000/WinPixDown*Tfactor*(0:256)/512;
    OrigPoints = length(vel);
    
    % plot everything
    fig1 = figure;
    set (fig1, 'Units', 'normalized', 'Position', [.1, .1, .8 .8] );
    Plotvel = subplot('Position',[.05 .7 .6 .25]);
    Plotflux = subplot('Position',[.05 .4 .6 .25]);
    Plotana = subplot('Position',[.05 .05 .6 .25]);
    Plothist = subplot('Position', [.7 .7 .25 .25]);
    Plotfreq = subplot('Position', [.7 .4 .25 .25]);
    Plotstim = subplot('Position', [.7 .05 .25 .25]);
    orient landscape;
    
    subplot(Plotvel);
    plot(time, vel, '-b.'); hold on;
    plot(time, smoothed, 'b')
    xlabel({['time (s)'], [Openfile]})
    ylabel('velocity (mm/s)')
    title(['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev)])
    hold on
    ys = ylim;
    
    subplot(Plotfreq)
    semilogy(f, Pyy(1:257))
    xlabel('Hz')
    hold on
    ys = ylim;
    
    subplot(Plothist);
    [n, xout] = hist(vel,25);
    barh(xout,n, 'b');
    ylim(ys)
    
    if use_flux % flux plot
        subplot(Plotflux);
        plot(time, flux1, 'r.', time, flux2, 'g.', time, flux3, 'b.');
        text(0.1, 0.9, ['avgs: ', num2str(avgflux1), '  ',num2str(avgflux2), '  ',num2str(avgflux3)])
        xlabel({['time (s)'], [Openfile]});
        ylabel('Flux (RBC/s)');
    end; % if plottin glfux
    
    
    if select_range
        
        good_range = 0;
        while not(good_range)
            subplot(Plotvel);
            plot(time, vel, '.');
            xlabel('time (s)')
            ylabel('velocity (mm/s)')
            title(['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev)])
            currentys = ylim;
            hold on
            
            
            subplot(Plotfreq)
            semilogy(f, Pyy(1:257))
            xlabel('Hz')
            hold on
            
                       
            if usestim~=0 %plot stimulus
                timestim = squeeze(Result(:,2))/1000;
                stim = squeeze(Result(:,usestim));
                figure(fig1); subplot(Plotstim)
                plot(timestim, stim, 'b'); hold on
                ylabel('stimulus')
                xlabel('time (ms)');
            end
            
            % ask user for input to select range of velocities
            set(gcf, 'Name', 'press a key to continue')
            w =0;
            while w == 0
                w = waitforbuttonpress;
            end

            limitok = 0;
            while not(limitok); % check that upper limit > lower limit
                prompt={'Upper Limit',...
                    'Lower Limit'}
                def={num2str(currentys(2)), num2str(currentys(1))};
                dlgTitle='Valid velocities for data';
                lineNo=1;
                answer=inputdlg(prompt,dlgTitle,lineNo,def,'on');
                UpLimit = str2double(answer(1));
                LowLimit = str2double(answer(2));
                if UpLimit > LowLimit
                    limitok = 1;
                else
                    beep
                end; %check limits ok 
            end; %while not(limitok)
            
            % Deal with Outliers
            outliers = (vel<=LowLimit) | (vel>=UpLimit);                        
            % remove outliers
            vel2 = vel;
            time2 = time;
            vel2(outliers) = [];
            time2(outliers) = [];
            if usestim~=0 %stimulus
                  stim2 = stim;
                  timestim2 = timestim;
                  stim2(outliers) = [];
                  timestim2(outliers) = [];
            end

            
            if use_flux 
                fluxtime12 = time2; fluxtime22 = time2; fluxtime32 = time2;
                flux12 = flux1;  flux22 = flux2;  flux32 = flux3;
                flux12(outliers) = []; flux22(outliers) = [];  flux32(outliers)= []; % remove velocity outliers 
                fluxtime12(flux12==10000) = []; flux12(flux12==10000) = []; % remove NaN fluxes
                fluxtime22(flux22==10000) = []; flux22(flux22==10000) = [];
                fluxtime32(flux32==10000) = []; flux32(flux32==10000) = [];
                fluxtime12(isnan(flux12)) = []; flux12(isnan(flux12)) = []; % remove NaN fluxes
                fluxtime22(isnan(flux22)) = []; flux22(isnan(flux22)) = [];
                fluxtime32(isnan(flux32)) = []; flux32(isnan(flux32)) = [];
                
                avgflux1 = mean(flux12); avgflux2= mean(flux22); avgflux3 = mean(flux32);
            end
            
            
            [NewPoints,p] = size(vel2);
            PercentValid = 100*NewPoints/OrigPoints;
            
            smoothed2 = convn(vel2, convbox, 'same');
            avg = mean(vel2);
            absavg = mean(abs(vel2));
            stdev = std(vel2);
            absstdev = std(abs(vel2));
            
            % Plot
            figure(fig1); subplot(Plotvel)
            plot(time2, vel2, '-ro'); hold on;
            plot(time2, smoothed2, 'r')
            xlabel({['time (s)'], [Openfile]})
            ylabel('velocity (mm/s)')
            title(['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev), ' % Valid:', num2str(PercentValid), '% ','Min: ', num2str(LowLimit), '  Max:', num2str(UpLimit)])
            hold off; ys = ylim;
            
            Y = fft(vel2, 512);
            Pyy = Y.* conj(Y) / 512;
            f = 1000/WinPixDown*Tfactor*(0:256)/512;
            subplot(Plotfreq)
            semilogy(f, Pyy(1:257), 'ro')
            xlabel('Hz')
            hold off
            
            figure(fig1)
            subplot(Plothist);
            [n2 xout2] = hist(vel2,25);
            barh(xout,n, 'b');
            hold on;
            barh(xout2, n2, 'r');
            hold off
            ylim(ys);
            title('PRESS A KEY WHEN DONE')
            
            
            if use_flux % flux plot
                figure(fig1)
                subplot(Plotflux);      
                plot(fluxtime12, flux12, 'r.', fluxtime22, flux22, 'g.', fluxtime32, flux32, 'b.');
                text(0.1, 0.9, ['avgs: ', num2str(avgflux1), '  ',num2str(avgflux2), '  ',num2str(avgflux3)])
                xlabel({['time (s)'], [Openfile]});
                ylabel('Flux (RBC/s)');
            end; % if plotting flux
            
            if usestim~=0 %plot stimulus for now just plot tract
                figure(fig1); subplot(Plotstim)
                plot(timestim2, stim2, 'b'); hold on
                ylabel('stimulus')
                xlabel('time (ms)');
                
                
                % % % % % % % % % % % % % % % % % % %
                % STIM-TRIGGERED AVERGAGE
                % GET STIM TIMES Tstim in seconds
                stimthresh = stim2 > threshforstim; %%%%%%%Need to get stim2 also
                Tstim = [];
                pointsleft = sum(stimthresh);
                count = 1;
                while pointsleft > 0
                    I = find(stimthresh);
                    Tstim(count) = time2(I(1));
                    toclear = I:I+(round((Tcycle-2) * Tfactor *1000)/WinPixDown);
                    stimthresh(toclear) = 0;
                    count = count +1;
                    pointsleft = sum(stimthresh)
                end
    
                Nstim = length(Tstim);
                Tstim
                
                % put appropriate numbers in bins
                timeCycle = -5:Tstep:Tcycle;
                velCycle = zeros(size(timeCycle));
                ecogCycle = zeros(size(timeCycle));
                nCycle = zeros(size(timeCycle));
                smoothedCycle = zeros(size(timeCycle));
                vels = zeros(Nstim, size(timeCycle));
                ns = zeros(Nstim, size(timeCycle));
                for k = 1:Nstim;
                    tstart = Tstim(k) - 5
                    tend = Tstim(k) +Tcycle
                    
                    InCycle = (time2>tstart) & (time2<tend);
                    testtime = time2(InCycle);
                    testvel = vel2(InCycle);
%                     testecog = ecog2(InCycle);
                    testsmoothed = smoothed2(InCycle);
                    
                    for j = 1:length(timeCycle);
                        InBin = ((testtime-Tstim(k)) >= (timeCycle(j)-Tbin)) & ((testtime-Tstim(k)) < timeCycle(j));
                        vels(k,j) = sum(testvel(InBin))/sum(InBin);
                    velCycle(j) = sum(testvel(InBin))/sum(InBin);
%                         ecogCycle(j) = sum(testecog(InBin))/sum(InBin);
                        smoothedCycle(j) = sum(testsmoothed(InBin))/sum(InBin);
                        nCycle(j) = sum(InBin);
                        ns(k,j) = sum(InBin);
                    end
                end
                
                vels(isnan(vels)) = 0;
                
                % Plot temp new fig
                subplot(Plotana);
                avgvels = sum(vels)./sum(ns);
                plot(timeCycle, avgvels); hold on;
                ylabel('stim trig avg')
                xlabel('cycle time')
                title(Openfile)
                
%                 subplot(Plotstim);
%                 hold off
%                 imagesc(vels); axis image;
%                 
%                 pause;
                
                save([pname  'stimtrig' files(i).name], 'vels', 'ns', 'timeCycle'); 
                
                % Cummulative for all data

                % % % % % % % % % % % % % % % % % % % 
                
          else
                %TEMP fornow plor flux vs.vel
%                 figure(fig1); subplot(Plotstim)
%                 flux1(outliers) = [];
%                 %                 plot(vel2, flux1, '.')
%                 %                 xlabel('velocity')
%                 %                 ylabel('flux1')
%                 vel2(vel2 == 0) = 0.00001;
%                 
%                 plot(time2, flux1/vel2)
%                 xlabel('time'); ylabel('density (RBC/mm)');
          end
          
          if useana~=0 % plot analog channel
                timeana = squeeze(Result(:,2))/1000;
                ana = squeeze(Result(:,useana));
                figure(fig1); subplot(Plotana)
                plot(timestim, ana, 'r'); hold on
                ylabel('analog')
                xlabel('time (ms)');
            end


            figure(fig1); subplot(Plotvel);
            ys = ylim;
            subplot(Plothist);
            ylim(ys);
            
            % Allow user to move stuff around
            set(gcf, 'Name', 'press a key to continue')
            w =0;
            while w == 0
                w = waitforbuttonpress;
            end

            button = questdlg('Use this range?',...
                'Continue Operation','Yes','No','No');
            if strcmp(button,'Yes')
                good_range = 1;
                
                %                                 %%%%%%%%%%%%%
                %                                 % find systolic, diastolic
                %                                 Tend = time2(end);
                %                                 
                %                                 use = zeros(size(time2));
                %                                 tstart = 0;
                %                                 
                %                                 systolic = [];
                %                                 diastolic = [];
                %                                 timesysdia = [];
                %                                 
                %                                 while tstart + Twin < Tend
                %                                         use = (time2>tstart) & (time2<(tstart+Twin));
                %                                         Max = max(vel2(use));
                %                                         Min = min(vel2(use));
                %                                         
                %                                         systolic = horzcat(systolic, Max);
                %                                         diastolic = horzcat(diastolic, Min);
                %                                         timesysdia = horzcat(timesysdia, tstart + Twin/2);
                %                                         
                %                                         
                %                                         tstart = tstart + WinPixDown/Tfactor/1000;
                %                                 end; % while tstart loop
                %                                 
                %                                 box = ones(1,  5*round(Twin*1000/WinPixDown))/5/round(Twin*1000/WinPixDown);
                %                                 
                %                                 systolic = convn(systolic, box, 'same');
                %                                 diastolic = convn(diastolic, box, 'same');
                %                                 
                %                                 systolic(1:5*round(Twin*1000/WinPixDown)) = [];
                %                                 systolic(length(systolic)-5*round(Twin*1000/WinPixDown):end) = [];
                %                                 diastolic(1:5*round(Twin*1000/WinPixDown)) = [];
                %                                 diastolic(length(diastolic)-5*round(Twin*1000/WinPixDown):end) = [];
                %                                 timesysdia(1:5*round(Twin*1000/WinPixDown)) = [];
                %                                 timesysdia(length(timesysdia)-5*round(Twin*1000/WinPixDown):end) = [];
                %                                 Avgsys = mean(systolic);
                %                                 Avgdia = mean(diastolic);
                %                                 %%%%%%%%%%%%%%% END PLOT SYS/DIASTOLIC
                %                                 
                %                                 subplot(Plotvel)
                %                                 hold on;
                %                                 plot(timesysdia(1:length(systolic)), systolic, 'k', timesysdia(1:length(diastolic)),diastolic, 'k');
                %                                 title({['Velocity: ', num2str(avg), 'mm/s  Std. Dev.: ', num2str(stdev), ' % Valid:', num2str(PercentValid), '%'],
                %                                         ['Min: ', num2str(LowLimit), '  Max:', num2str(UpLimit), ' Sys Vel:', num2str(Avgsys), ' Dia Vel:', num2str(Avgdia)] }, 'FontSize', 8)
                %                                 
                
%                 close fig1
                
            end               
            
        end
    end
    
    Savefile = strrep(Openfile, '.mat', 'ExcludedPts.mat')
    save(Savefile, 'time2', 'vel2');
    
    if use_many
        Names = strvcat(Names, files(i).name);
        Velocities = strvcat(Velocities, num2str(avg));
        StdDev = strvcat(StdDev, num2str(stdev));
        AbsVelocities = strvcat(AbsVelocities, num2str(absavg));
        AbsStdDev = strvcat(AbsStdDev, num2str(absstdev));
        if select_range
            MinRanges = strvcat(MinRanges, num2str(LowLimit));
            MaxRanges = strvcat(MaxRanges, num2str(UpLimit));
            PerValids = strvcat(PerValids, num2str(PercentValid));
%             AvgSyss = strvcat(AvgSyss, num2str(Avgsys));
%             AvgDias = strvcat(AvgDias, num2str(Avgdia));
            %                 close (fig1);
        end
        if use_flux
            AvgFlux1s  = strvcat(AvgFlux1s, num2str(avgflux1));
            AvgFlux2s  = strvcat(AvgFlux2s, num2str(avgflux2));
            AvgFlux3s  = strvcat(AvgFlux3s, num2str(avgflux3));
        end
        
        % save cummulative data everytime
        commas = [];
        for j = 1:i+1
            commas(j,1) = ',';
        end
        if select_range
            tosave = horzcat(char(Names), commas, char(Velocities), commas, char(StdDev), commas,char(AbsVelocities), commas, char(AbsStdDev),commas, MinRanges, commas, MaxRanges,commas, PerValids); %,commas, AvgSyss, commas, AvgDias); TEMP
        else
            tosave = horzcat(char(Names), commas, char(Velocities), commas, char(StdDev), commas,char(AbsVelocities), commas,char(AbsStdDev));
        end
        if use_flux
            tosave = horzcat(tosave, commas, char(AvgFlux1s),commas, char(AvgFlux2s),commas, char(AvgFlux3s));
        end

        diary(Datafiletemp)
        tosave(i+1,:)
        diary off;
        
        Savefile = strrep(Openfile, '.mat', 'ExcludedPts.mat')
        save(Savefile, 'time2', 'vel2');

    end; % if use many
    
   close(fig1)
end; % loop through movies



if use_many % write data

    commas = [];
    for i = 1:nmovies+1
        commas(i,1) = ',';
    end
    if select_range
        tosave = horzcat(char(Names), commas, char(Velocities), commas, char(StdDev), commas,char(AbsVelocities), commas,char(AbsStdDev), commas, MinRanges, commas, MaxRanges,commas, PerValids); %,commas, AvgSyss, commas, AvgDias); TEMP
    else
        tosave = horzcat(char(Names), commas, char(Velocities), commas, char(StdDev), commas,char(AbsVelocities), commas,char(AbsStdDev));
    end
    if use_flux
        tosave = horzcat(tosave, commas, char(AvgFlux1s),commas, char(AvgFlux2s),commas, char(AvgFlux3s));
    end

    diary(Datafile)
    tosave
    diary off
end
