function next_mcs=predict(SNR_history,threshold,upper_threshold,lower_threshold,current_mcs,mode)

    next_mcs=current_mcs;

    window_size=3;

    history_weight=[0.1^2,0.1,1-0.1-0.1^2];

    threshold=[-inf,threshold,inf];

    if mode=="sliding_and_change"

        if (length(SNR_history)>window_size)
            SNR_predicting=sum(history_weight.*SNR_history([end-window_size+1:1:end]));
            increaseMCS = (SNR_predicting > upper_threshold((current_mcs==0)+current_mcs));
            decreaseMCS = (SNR_predicting <= lower_threshold((current_mcs==0)+current_mcs));
            next_mcs = current_mcs+increaseMCS-decreaseMCS;

        else
            SNR_predicting=SNR_history(end);
            increaseMCS = (SNR_predicting > upper_threshold((current_mcs==0)+current_mcs));
            decreaseMCS = (SNR_predicting <= lower_threshold((current_mcs==0)+current_mcs));
            next_mcs = current_mcs+increaseMCS-decreaseMCS; 
        end


    elseif mode=="sliding_and_choose"

        if (length(SNR_history)>window_size)
            % SNR_predicting=sum(history_weight.*SNR_history([end-window_size+1:1:end]));
            % SNR_predicting=SNR_history(end)+(SNR_history(end)-SNR_history(end-1));
            SNR_predicting=SNR_history(end);
            for i=1:1:length(threshold)-1
                
                % if SNR_predicting>=lower_threshold(i) && SNR_predicting<upper_threshold(i)

                if SNR_predicting>=threshold(i) && SNR_predicting<threshold(i+1)
                    
                    next_mcs=i;

                    break;

                end
            
            end

        else

            SNR_predicting=SNR_history(end);

            for i=1:1:length(threshold)-1
                if SNR_predicting>=threshold(i) && SNR_predicting<threshold(i+1)
                    
                    next_mcs=i;

                    break;

                end
            end
        end

    % elseif mode=="ARIMA"

    %     if (length(SNR_history)>50)
    %         SNR_predicting= ARIMA_predict(SNR_history',1,3,3,'off');

    %         for i=1:1:length(threshold)-1
                
    %             % if SNR_predicting>=lower_threshold(i) && SNR_predicting<upper_threshold(i)

    %             if SNR_predicting>=threshold(i) && SNR_predicting<threshold(i+1)
                    
    %                 next_mcs=i;

    %                 break;

    %             end
            
    %         end
    %     else
    %         SNR_predicting=SNR_history(end);

    %         for i=1:1:length(threshold)-1
    %             if SNR_predicting>=threshold(i) && SNR_predicting<threshold(i+1)
                    
    %                 next_mcs=i;

    %                 break;

    %             end
    %         end
    %     end



    end





    % if (next_mcs<1)
    %     next_mcs=1;
    % elseif (next_mcs>10)
    %     next_mcs=10;
    % end

           
end