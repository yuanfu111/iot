
function results=Rate_Control_Algorithm_1(tgacChannel,cfgVHT,numPackets,baseSNR,maxJump,snrWalk,meanSNR,amplitude,sr,predict_mode)

    % Set random stream for repeatability of results


    s = rng(21); % It will be restored later as rng(s)
    %% Rate Control Algorithm Parameters
    %
    % From here - These are input parameters for an example rate control algorithm
    % [PLEASE REPLACE IT WITH YOUR OWN SET OF PARAMETERS/CHANGE BASED ON YOUR DESIGN]
    %
    rcaAttack = 1;  % Control the sensitivity when MCS is increasing
    rcaRelease = 0; % Control the sensitivity when MCS is decreasing
    threshold = [11 14 19 20 25 28 30 31 35];
    % threshold=[10:3:40];
    snrUp = [threshold inf] + rcaAttack;
    snrDown = [-inf threshold] - rcaRelease;
    % expected_speed=[]
    snrInd = cfgVHT.MCS; % Store the start MCS value
    %
    % [PLEASE REPLACE IT WITH YOUR OWN SET OF PARAMETERS/CHANGE BASED ON YOUR DESIGN]
    SNR_history=[];

    % long-term channel condition evaluate
    channel_windows_size=30;
    current_channel_ber=0;
    previous_channel_ber=current_channel_ber;
    current_channel_per=0;
    current_table=threshold;
    current_upper_bound=snrUp;
    current_lower_bound=snrDown;
    threshold_changing_boundary=4;
    % MCS_in_Window=zeros(1,channel_windows_size);

    % short-term signal trend evaluate
    signal_windows_size=5;
    current_signal_ber=0;
    previous_signal_ber=current_signal_ber;
    current_signal_per=0;
    % current_signal_difference=0; %signal difference
    % previous_signal_per=current_signal_per;
    current_bias=0;
    bias_boundary=2;


    % To here - These are input parameters for an example rate control algorithm

    %% Define simulation variables

    snrMeasured = zeros(1,numPackets);
    MCS = zeros(1,numPackets);
    ber = zeros(1,numPackets);
    packetLength = zeros(1,numPackets);

    %% Processing Chain

    for numPkt = 1:numPackets 

        snrWalk = 0.9 * snrWalk + 0.1 * baseSNR(numPkt) + rand(1) * maxJump * 2 - maxJump; % DO NOT CHANGE THIS % Generate SNR value per packet using random walk algorithm biased towards the mean SNR
        
        txPSDU = randi([0,1],8 * cfgVHT.PSDULength, 1, 'int8');        % DO NOT CHANGE THIS % randomly generate bits of cfgVHT.PSDULength length
        txWave = wlanWaveformGenerator(txPSDU,cfgVHT,'IdleTime',5e-4); % DO NOT CHANGE THIS % randomly generate a single packet waveform
        
        y = processPacket(txWave,snrWalk,tgacChannel,cfgVHT);          % DO NOT CHANGE THIS % Receive processing, including SNR estimation
        
        % Store estimated SNR value for each packet
        if isempty(y.EstimatedSNR) 
            snrMeasured(1,numPkt) = NaN;
%             SNR_history=[SNR_history,0];
        else
            snrMeasured(1,numPkt) = y.EstimatedSNR;
%             SNR_history=[SNR_history,y.EstimatedSNR];
        end
        
        % Determine the length of the packet in seconds including idle time
        packetLength(numPkt) = y.RxWaveformLength/sr;
        
        % Calculate packet error rate (PER)
        if isempty(y.RxPSDU)
            % Set the PER of an undetected packet to NaN
            ber(numPkt) = NaN;
        else
            [~,ber(numPkt)] = biterr(y.RxPSDU,txPSDU);
        end
        
        % From here - This is an example rate control algorithm
        % [PLEASE REPLACE IT WITH YOUR OWN ALGORITHM]
        
        % Store estimated SNR value for each packet
        if isempty(y.EstimatedSNR) 
            SNR_history=[SNR_history,0];
        else
            SNR_history=[SNR_history,y.EstimatedSNR];
        end

        if (ber(numPkt)<=1)
            current_channel_ber=current_channel_ber+ber(numPkt)/channel_windows_size;
            current_signal_ber=current_signal_ber+ber(numPkt)/signal_windows_size;
            if (ber(numPkt)>0)
                current_signal_per=current_signal_per+1;
                current_channel_per=current_channel_per+1;
            end
            
        else 
            current_channel_ber=current_channel_ber+1/channel_windows_size;
            current_signal_ber=current_signal_ber+1/signal_windows_size;
            current_signal_per=current_signal_per+1;
            current_channel_per=current_channel_per+1;
        end

        % update bias

        if (mod(numPkt,signal_windows_size)==0)

            current_signal_ber=0.1*previous_signal_ber+0.9*current_signal_ber;

            previous_signal_ber=current_signal_ber;
            
                if (current_signal_ber==0.000) || (current_signal_per<=0.1*signal_windows_size) 

                    if (current_bias<0)

                        current_bias=current_bias+1;

                    elseif (current_bias<bias_boundary) && (current_table(1)-threshold(1)>0)

                        current_bias=current_bias+1;

                    elseif (current_bias<bias_boundary-1) && (current_table(1)-threshold(1)==0)
                        current_bias=current_bias+1;
                    end

                elseif (current_signal_ber>0.01) || (current_signal_per>0.2*signal_windows_size)

                    if (current_bias>-bias_boundary) 

                        current_bias=current_bias-1;

                    end

                end
            % end

            current_signal_ber=0;
            current_signal_per=0;
        end

        %  update table/threshold
        if (mod(numPkt,channel_windows_size)==0)
            current_channel_ber=0.1*previous_channel_ber+0.9*current_channel_ber;
            current_channel_ber;
            current_channel_per;
            previous_channel_ber=current_channel_ber;
            [current_table,current_upper_bound,current_lower_bound]=compensate(threshold,threshold_changing_boundary,current_table,current_channel_ber,current_channel_per, channel_windows_size);
                current_table;
            % end
            current_channel_ber=0;
            current_channel_per=0;
            % current_bias=0;
            % MCS_in_Window=zeros(1,channel_windows_size);
        end


        MCS(numPkt) = cfgVHT.MCS;
        current_mcs=cfgVHT.MCS;
        next_mcs=predict(SNR_history,threshold,current_upper_bound,current_lower_bound,current_mcs,predict_mode);


        if ((next_mcs+current_bias-1<=9) && (next_mcs+current_bias-1>=0))

            next_mcs=next_mcs+current_bias;
        
        end

        cfgVHT.MCS=next_mcs-1;


        % [PLEASE REPLACE IT WITH YOUR OWN ALGORITHM]
        % To here - This is an example rate control algorithm
    end
    % ber
    %% Display and Plot Simulation Results

    % Display and plot simulation results
    ov_dr = num2str(8*cfgVHT.APEPLength*(numPackets-numel(find(ber)))/sum(packetLength)/1e6);
    disp(['Overall data rate: ' ov_dr ' Mbps']);
    ov_per = num2str(numel(find(ber))/numPackets);
    disp(['Overall packet error rate: ' ov_per]);

    plotResults(ber,packetLength,snrMeasured,MCS,cfgVHT);
    current_upper_bound
    current_lower_bound
    % SNR_history;

    % Restore default stream
    rng(s);


    results = struct( ...
    'meanSNR',          meanSNR, ...
    'amplitude',        amplitude, ...
    'maxJump',          maxJump, ...
    'distance',         tgacChannel.TransmitReceiveDistance, ...
    'numPackets',       numPackets, ....
    'mode',             predict_mode, ...
    'window_size',      channel_windows_size, ...
    'Data_rate',        ov_dr, ...
    'PER',              ov_per);

end



% NOTE: TRY NOT TO CHANGE THIS FUNCTION - let us make the packet processing the same for everyone
function Y = processPacket(txWave,snrWalk,tgacChannel,cfgVHT)
    % Pass the transmitted waveform through the channel, perform
    % receiver processing, and SNR estimation.
    
    chanBW = cfgVHT.ChannelBandwidth; % Channel bandwidth
    % Set the following parameters to empty for an undetected packet
    estimatedSNR = [];
    eqDataSym = [];
    noiseVarVHT = [];
    rxPSDU = [];
    
    ofdmInfo = wlanVHTOFDMInfo('VHT-Data',cfgVHT); % Get the OFDM info
    rxWave = tgacChannel(txWave); % Pass the waveform through the fading channel model
   
    awgnChannel = comm.AWGNChannel; % Create an instance of the AWGN channel for each transmitted packet
    awgnChannel.NoiseMethod = 'Signal to noise ratio (SNR)';
    awgnChannel.SignalPower = 1/tgacChannel.NumReceiveAntennas; % Normalization
    awgnChannel.SNR = snrWalk - 10 * log10(ofdmInfo.FFTLength/ofdmInfo.NumTones); % Account for energy in nulls
    
    rxWave = awgnChannel(rxWave); % Add noise
    rxWaveformLength = size(rxWave,1); % Length of the received waveform
    
    % Recover packet
    ind = wlanFieldIndices(cfgVHT); % Get field indices
    pktOffset = wlanPacketDetect(rxWave,chanBW); % Detect packet
    
    if ~isempty(pktOffset) % If packet detected
        LLTFSearchBuffer = rxWave(pktOffset+(ind.LSTF(1):ind.LSIG(2)),:); % Extract the L-LTF field for fine timing synchronization
        finePktOffset = wlanSymbolTimingEstimate(LLTFSearchBuffer,chanBW); % Start index of L-LTF field
        pktOffset = pktOffset + finePktOffset; % Determine final packet offset
        
        if pktOffset < 15 % If synchronization successful
            % Extract VHT-LTF samples from the waveform, demodulate and perform channel estimation
            VHTLTF = rxWave(pktOffset + (ind.VHTLTF(1):ind.VHTLTF(2)),:);
            demodVHTLTF = wlanVHTLTFDemodulate(VHTLTF,cfgVHT);
            chanEstVHTLTF = wlanVHTLTFChannelEstimate(demodVHTLTF,cfgVHT);
            chanEstSSPilots = vhtSingleStreamChannelEstimate(demodVHTLTF,cfgVHT); % Get single stream channel estimate
            vhtdata = rxWave(pktOffset+(ind.VHTData(1):ind.VHTData(2)),:); % Extract VHT data field
            noiseVarVHT = vhtNoiseEstimate(vhtdata,chanEstSSPilots,cfgVHT); % Estimate the noise power in VHT data field
            
            % Recover equalized symbols at data carrying subcarriers using channel estimates from VHT-LTF
            [rxPSDU,~,eqDataSym] = wlanVHTDataRecover(vhtdata,chanEstVHTLTF,noiseVarVHT,cfgVHT);
            
            % SNR estimation per receive antenna
            powVHTLTF = mean(VHTLTF.*conj(VHTLTF));
            estSigPower = powVHTLTF - noiseVarVHT;
            estimatedSNR = 10*log10(mean(estSigPower./noiseVarVHT));
        end
    end
    
    % Set output
    Y = struct( ...
        'RxPSDU',           rxPSDU, ...
        'EqDataSym',        eqDataSym, ...
        'RxWaveformLength', rxWaveformLength, ...
        'NoiseVar',         noiseVarVHT, ...
        'EstimatedSNR',     estimatedSNR);
    
end

function plotResults(ber,packetLength,snrMeasured,MCS,cfgVHT)
    % Visualize simulation results

    figure('Outerposition',[50 50 900 700])
    subplot(4,1,1);
    plot(MCS);
    xlabel('Packet Number')
    ylabel('MCS')
    title('MCS selected for transmission')

    subplot(4,1,2);
    plot(snrMeasured);
    xlabel('Packet Number')
    ylabel('SNR')
    title('Estimated SNR')

    subplot(4,1,3);
    plot(find(ber==0),ber(ber==0),'x') 
    hold on; stem(find(ber>0),ber(ber>0),'or') 
    if any(ber)
        legend('Successful decode','Unsuccessful decode') 
    else
        legend('Successful decode') 
    end
    xlabel('Packet Number')
    ylabel('BER')
    title('Instantaneous bit error rate per packet')

    subplot(4,1,4);
    windowLength = 3; % Length of the averaging window
    movDataRate = movsum(8*cfgVHT.APEPLength.*(ber==0),windowLength)./movsum(packetLength,windowLength)/1e6;
    plot(movDataRate)
    xlabel('Packet Number')
    ylabel('Mbps')
    title(sprintf('Throughput over last %d packets',windowLength))
    
end

