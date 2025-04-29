% UV Method (MODI Method) for optimisation of Transportation Problem

% Taking the inputs, assuming the IBFS and updated cost matrix will be
% given
allocations = input("Enter the final allocations matrix obtained as IBFS: ");
cost = input("Enter the updated cost matrix: ");
disp("These are current allocations: ");
disp(allocations);
[num_sources, num_destinations] = size(cost);

while true
    %first we will initialise the u_i and v_j, ofc as 0 as of now
    % and then will put through the loop to update
    u = zeros(num_sources, 1);
    v = zeros(1, num_destinations);
    u_filled = false(num_sources, 1);
    v_filled = false(1, num_destinations);

    u(1) = 0;    % we are letting the u(1) = 0 as we need to assume one of these variables
    % as there are m+n-1 equations ( corresponding to the basic vars) but
    % m+n variables (u_i and v_j)
    u_filled(1) = true;
    
    % now the loop calculates the u_i and v_j
    % using the 
    while ~all(u_filled) || ~all(v_filled)
        for i = 1:num_sources
            for j = 1:num_destinations
                if allocations(i,j) > 0
                    if u_filled(i) && ~v_filled(j)
                        v(j) = cost(i,j) - u(i);
                        v_filled(j) = true;
                    elseif v_filled(j) && ~u_filled(i)
                        u(i) = cost(i,j) - v(j);
                        u_filled(i) = true;
                    end
                end
            end
        end
    end

    % now we will calculate the u_i+v_j-c_i_j for non basic cells 
    % which will help in making the decision that we have reached the
    % optimal solution or we need to process further
    % and which cell to process
    opportunity_cost = zeros(num_sources, num_destinations);
    for i = 1:num_sources
        for j = 1:num_destinations
            opportunity_cost(i,j) = u(i) + v(j) - cost(i,j);
        end
    end

    % ofc entering variable would be the one with most positive from non
    % basic cells, for basic it would be 0 only
    [max_val, idx] = max(opportunity_cost(:));
    [enter_row, enter_col] = ind2sub(size(opportunity_cost), idx);

    % if max is 0 , ie for all non basic cells its -ve only 
    % thus optimal assignment is reached and we can take an exist here
    if max_val <= 0
        disp('Optimal solution found.');
        break;
    end

    % now we need to process further , as it didnt take an exit above
    % we have the entering cell , we need outgoing cell and adjust the
    % assignments as well
    % this can only handle rectangle loops
    loop_found = false;
    for i = 1:num_sources
        for j = 1:num_destinations
            % looking for the loop here wrt the entering cell
            if allocations(i,enter_col) > 0 && allocations(enter_row,j) > 0 && allocations(i,j) > 0
                exit_row = i;
                exit_col = j;
                loop_found = true;
                break;
            end
        end
        if loop_found
            break;
        end
    end
    % random error handling
    if ~loop_found
        error('No adjustment loop found. IBFS may not be correct.');
    end

    % now once we have the loop corners then we assign '+' to entering and
    % opposite to it and rest two '-' 
    % and we get the theta from the '-' assinged corners and update
    % accoringly
    if allocations(exit_row, enter_col) < allocations(enter_row, exit_col)
        theta = allocations(exit_row, enter_col);
        
        allocations(exit_row, enter_col) = 0;
        allocations(enter_row, exit_col) = allocations(enter_row, exit_col) - theta;
        allocations(enter_row, enter_col) = theta;
        allocations(exit_row, exit_col) = allocations(exit_row, exit_col) + theta;
    else
        theta = allocations(enter_row, exit_col);
        
        allocations(enter_row, exit_col) = 0;
        allocations(exit_row, enter_col) = allocations(exit_row, enter_col) - theta;
        allocations(enter_row, enter_col) = theta;
        allocations(exit_row, exit_col) = allocations(exit_row, exit_col) + theta;
    end
    disp("Updated alloctions: ");
    disp(allocations);
end

% now once we are out of loop we have the optimal assignment and we can get
% the final cost/value of the optimal assignment
final_cost = 0;
for i = 1:num_sources
    for j = 1:num_destinations
        final_cost = final_cost + allocations(i,j) * cost(i,j);
    end
end
disp('Final Allocation Matrix:');
disp(allocations);
disp(['Min Transportation Cost: ', num2str(final_cost)]);
