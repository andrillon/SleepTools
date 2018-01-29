% For example
% [upStateTimes, downStateTimes] = detect_SW_negative_withStage_segments_dealWithSAW_scalp_new('REC1M_EEG_BP.mat', 1000, 2000, 'stagingRevised.mat', 5);
function [allWaves, slowWaves] = SWsDetectionAlgorithm_forSPM(D, SleepScoring_filename, paramDetection)
% (EEG_filename, original_fs, desiredPercentOfSlowWaves, stagingFile, SAW_threshold, onlyThisStage)
allWaves = [];
slowWaves=[];
SR=D.fsample;

if nargin<3
    paramDetection=[];
    fprintf('\n *** Parameters set to default ***\n')
end
if isfield(paramDetection,'SWband'),              SWband=paramDetection.SWband;                           else SWband=[0.2 3 3]; end
if isfield(paramDetection,'ChannelSelection'),    ChannelSelection=paramDetection.ChannelSelection;       else ChannelSelection=D.meegchannels; end
if isfield(paramDetection,'SWslope'),             SWslope=paramDetection.SWslope;                         else SWslope=[0.25 1]; end
if isfield(paramDetection,'P2Pamp'),              P2Pamp=paramDetection.P2Pamp;                           else P2Pamp=75; end

%% Filter EEG
if length(SWband)==3
    S=[];
    S.D=D;
    S.filter.type='but';
    S.filter.order=SWband(3);
    S.filter.band='bandpass';
    S.filter.PHz=SWband;
    S.filter.dir='twopass';
    S.save=0;
    Dfilt=spm_eeg_filter(S);
else
    S=[];
    S.D=D;
    S.filter.type='but';
    S.filter.order=SWband(2);
    S.filter.band='lowpass';
    S.filter.PHz=SWband(1);
    S.filter.dir='twopass';
    S.save=0;
    Dfilt=spm_eeg_filter(S);
end
fsample_ori=D.fsample;

%% Down sample to 100Hz
fsample_new=100;
S=[];
S.D=Dfilt;
S.fsample_new=100;
S.save=0;
Dres=spm_eeg_downsample(S);
SR=Dres.fsample;

countChan=0;
for nChan=ChannelSelection
    countChan=countChan+1;
    fprintf('... detecting SWs on %s (%g/%g)',Dres.chanlabels{nChan},countChan,length(ChannelSelection))
    EEGdata_BP_ss=squeeze(Dres(nChan,:,:));
    %     EEGdata_BP=bandpass(EEGdata, SR, SWband(1), SWband(2), SWband(3));
    
    % Include sleep stagin
    if isempty(SleepScoring_filename)
        SleepStages_ts=nan(1,length(EEGdata_BP_ss));
    else
        [SleepStages_ts arousal_sampled marousal_sampled]=format_sleepScoring(Dres,SleepScoring_filename);
    end
    
    %%%%%% Brady's code from here on
    % EEG is the data vector sampled at 100Hz
    pos_index=zeros(length(EEGdata_BP_ss(1, :)),1);
    pos_index(find(EEGdata_BP_ss(1, :)>0))=1; %index of all positive points for EEG
    difference=diff(pos_index); poscross=find(difference==1) ; negcross=find(difference==-1); %find neg ZX and pos ZX
    EEGder=meanfilt(diff(EEGdata_BP_ss(1, :)),5);
    
    pos_index=zeros(length(EEGder),1);
    pos_index(find(EEGder>0.1))=1; %index of all positive points above minimum threshold
    difference=diff(pos_index);
    peaks=find(difference==-1)+1; troughs=find(difference==1)+1; %find pos ZX and neg ZX of the derivative (the peaks & troughs)
    peaks(EEGdata_BP_ss(1, peaks)<0)=[];troughs(EEGdata_BP_ss(1, troughs)>0)=[]; % rejects peaks below zero and troughs above zero
    
    if isempty(negcross) || isempty(poscross)
        allWaves{nChan}  =[];
        slowWaves{nChan} =[];
        fprintf('... no wave found\n')
        continue;
    end
    if negcross(1)<poscross(1);start=1;else start=2;end %makes negcross and poscross same size to start
    if start==2;poscross(1)=[];end
    
    lastpk=NaN; %way to look at Peak to Peak parameters if needed
    
    waves = zeros(length(negcross)-start, 23);
    
    for wndx=start:length(negcross)-1
        
        wavest=negcross(wndx);   %only used for neg/pos peaks
        wavend=negcross(wndx+1); %only used for neg/pos peaks
        mxdn=abs(min(meanfilt(diff(EEGdata_BP_ss(1, wavest:poscross(wndx))),5)))*SR;     % matrix (27) determines instantaneous positive 1st segement slope on smoothed signal, (name not representative)
        mxup=max(meanfilt(diff(EEGdata_BP_ss(1, wavest:poscross(wndx))),5))*SR;     % matrix (28) determines maximal negative slope for 2nd segement (name not representative)
        negpeaks=troughs(troughs>wavest&troughs<wavend);
        
        % In case a peak is not detected for this wave (happens rarely)
        if (size(negpeaks,1) == 0)
            waves(wndx, :) = NaN; %uvValueLine;
            continue;
        end
        
        pospeaks=peaks(peaks>wavest&peaks<=wavend);
        if isempty(pospeaks);pospeaks=wavend; end %if negpeaks is empty set negpeak to pos ZX
        period=wavend-wavest; %matrix(11) /SR
        poszx=poscross(wndx); %matrix(10)
        b=min(EEGdata_BP_ss(1, negpeaks)); % matrix (12) most pos peak /abs for matrix
        if b>0;b=b(1);end;
        bx=negpeaks(EEGdata_BP_ss(1, negpeaks)==b); %matrix (13) max pos peak location in entire night
        c=max(EEGdata_BP_ss(1, pospeaks)); % matrix (14) most neg peak
        if c>0;c=c(1);end;
        cx=pospeaks(EEGdata_BP_ss(1, pospeaks)==c); %matrix (15) max neg peak location in entire night
        maxb2c=c-b; % %matrix (16) max peak to peak amp
        nump=length(negpeaks); %matrix(24) now number of positive peaks
        n1=abs(EEGdata_BP_ss(1, negpeaks(1))); %matrix(17) 1st pos peak amp
        n1x=negpeaks(1); %matrix(18) 1st pos peak location
        nEnd=abs(EEGdata_BP_ss(1, negpeaks(end))); %matrix(19) last pos peak amp
        nEndx=negpeaks(end);%matrix(20) last pos peak location
        p1=EEGdata_BP_ss(1, pospeaks(1)); %matrix(21) 1st neg peak amp
        p1x=pospeaks(1); %matrix(22) 1st pos peak location
        meanAmp=abs(mean(EEGdata_BP_ss(1, negpeaks))); %matrix(23)
        nperiod=poszx-wavest; %matrix (25)neghalfwave period
        mdpt=wavest+ceil(nperiod/2); %matrix(9)
        p2p=(cx-lastpk)/SR; %matrix(26) 1st peak to last peak period
        lastpk=cx;
        
        thisStage=SleepStages_ts(poszx);
        % Result Matrix
        %1:  wave beginning (sample)
        %2:  wave end (sample)
        %3:  wave middle point (sample)
        %4:  wave negative half-way (sample)
        %5:  period in seconds
        %6:  positive amplitude peak
        %7:  positive amplitude peak position (sample)
        %8:  negative amplitude peak
        %9:  negative amplitude peak position (sample)
        %10: peak-to-peak amplitude
        %11: 1st pos peak amplitude
        %12: 1st pos peak amplitude position (sample)
        %13: Last pos peak amplitude
        %14: Last pos peak amplitude position (sample)
        %15: 1st neg peak amplitude
        %16: 1st neg peak amplitude position (sample)
        %17: mean wave amplitude
        %18: number of positive peaks
        %19: wave negative half-way period
        %20: 1st peak to last peak period
        %21: determines instantaneous positive 1st segement slope on smoothed signal
        %22: determines maximal negative slope for 2nd segement
        %23: stage (if scored data)
        waves(wndx, :) = [wavest wavend mdpt poszx period/SR abs(b) bx c cx maxb2c n1 n1x nEnd nEndx p1 p1x meanAmp nump nperiod/SR p2p mxdn mxup thisStage];
        waves(wndx, [1 2 3 4 7 9 12 14 16])=waves(wndx, [1 2 3 4 7 9 12 14 16])*(fsample_ori/fsample_new);

    end 
    allWaves{nChan} = waves;
 % of loop through segments

slowWaves{nChan} = waves((waves(:, 19)<SWslope(2)) & (waves(:, 19)>SWslope(1)) & (waves(:,10)>P2Pamp), :); % choose slow waves based on their period
fprintf('... %g SWs detected',size(slowWaves{nChan},1))
fprintf('\n')
end
end

function smdata=meanfilt(data,smwin)

if size(data,2)>smwin
    temp=nan(smwin,size(data,2)); % uses a 5 sample moving window to smooth derivative
    for n=1:smwin
        temp(n,:)=[nan(1,n-1) data(1,1:length(data)-n+1)];
    end
    smdata=mean(temp);
else
    smdata=data;
end
end