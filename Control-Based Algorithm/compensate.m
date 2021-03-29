function [new_table,new_upper_bound,new_lower_bound]=compensate(threshold,threshold_changing_boundary,current_table, channel_ber,channel_per,window_size)

    new_table=current_table;
    

    if (channel_ber<0.01) || (channel_per<0.1*window_size)
        % channel situation is better than expectation
        if threshold(1)-current_table(1)<threshold_changing_boundary
            new_table=new_table-1.0*ones(1,length(new_table));
        end

    elseif ((channel_ber>0.1) || (channel_per>0.2*window_size))
        % channel worse than expectation

        % length(find(MCS_in_Window==0))
        if current_table(1)-threshold(1)<threshold_changing_boundary
            new_table=new_table+1.0*ones(1,length(new_table));  
        end
    
    end

    new_upper_bound = [new_table inf] + 1;
    new_lower_bound = [-inf new_table] - 0;

end