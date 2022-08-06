close all;
clc;
clear;
% start reading input file
file_name = get_input_file();
fileID = fopen(file_name,'r');
buffer = fscanf(fileID,'%f');
fclose(fileID);
index = 5;
distance_mode = buffer(index, 1); % 1 is euclidean distance, 2 is manhattan distance
index = index + 1;
city_num = buffer(index,1);
index = index + 1;
map = zeros(city_num,2);
for i = 1 : city_num
    map(i,:) = [buffer(index,1),buffer(index+1,1)];
    index = index + 2;
end
width = [min(map,[],  1) - 1; max(map, [], 1) + 1];
route = asymmetric_traveling_salesman(map, distance_mode);
plot_final_graph(route, width);

function route = asymmetric_traveling_salesman(map, distance_mode)
    map_length = size(map,1);
    data = zeros(map_length);
    for i = 1 : map_length
        data(:, i) = calculate_distance(map, map(i, :), distance_mode);
    end
    max_weight = max(data, [], 'all');
    for i = 1 : map_length
        data(i, i) = Inf;
    end
    infinite_factor = 100;
    error_rate = 0.01;
    break_count = 1000;
    A = max_weight;
    init_t = max_weight;
    min_t = (error_rate/(map_length * (log2(map_length) * 3 / 2)));
    beta = 1.0 / ( max_weight * infinite_factor);
    same_energy_count = 0;
    current_t=init_t;
    spin = false(map_length);
    first_sum_v = map_length; 
    first_sum_j = zeros(map_length,1);
    second_sum_j = map_length;
    second_sum_v = zeros(map_length,1);
    third_sum = 0;
    first_sum_j_next = zeros(map_length,1);
    second_sum_v_next = zeros(map_length,1);
    total_energy = 0;
    last_total_energy = Inf;
    while(current_t>min_t)            
        for v = 1 : map_length
            for j = 1 : map_length
                if(spin(v,j))% 1 to 0
                    first_sum_j_next(v,1) = first_sum_j(v,1) - 1;
                    second_sum_v_next(j,1) = second_sum_v(j,1) - 1; 
                    third_difference = asymmetric_traveling_salesman_calculate_third_difference(spin, data, v, j, false) ;                                           
                else % 0 to 1
                    first_sum_j_next(v,1) = first_sum_j(v,1) + 1;
                    second_sum_v_next(j,1) = second_sum_v(j,1) + 1;
                    third_difference = asymmetric_traveling_salesman_calculate_third_difference(spin, data, v, j, true) ;     
                end
                first_sum_v_next = first_sum_v - mpower((1-first_sum_j(v,1)),2) + mpower((1-first_sum_j_next(v,1)),2);
                second_sum_j_next = second_sum_j - mpower((1-second_sum_v(j,1)),2) + mpower((1-second_sum_v_next(j,1)),2);
                energy_difference = (A * (first_sum_v_next - first_sum_v + second_sum_j_next - second_sum_j)) + third_difference;
                if (energy_difference < 0)
                    spin(v,j) = (~spin(v,j));
                    first_sum_v = first_sum_v_next; 
                    first_sum_j(v,1) = first_sum_j_next(v,1);
                    second_sum_j = second_sum_j_next;
                    second_sum_v(j,1) = second_sum_v_next(j,1);
                    third_sum = third_sum + third_difference;
                    total_energy = first_sum_v + second_sum_j +third_sum;                       
                else
                    prob = exp(-energy_difference/current_t);
                    if (prob>rand)
                        spin(v,j) = (~spin(v,j));
                        first_sum_v = first_sum_v_next; 
                        first_sum_j(v,1) = first_sum_j_next(v,1);
                        second_sum_j = second_sum_j_next;
                        second_sum_v(j,1) = second_sum_v_next(j,1);
                        third_sum = third_sum + third_difference;
                        total_energy = first_sum_v + second_sum_j +third_sum; 
                    end                        
                end      
            end                
        end
        fprintf("first sum  %d second sum %d third sum %d current t %f total energy %d\n",first_sum_v, second_sum_j, third_sum, current_t, total_energy);
        current_t = current_t / (1 + beta * current_t); 
        if (last_total_energy == total_energy)
            if (same_energy_count == break_count)
                break;
            end            
            same_energy_count = same_energy_count + 1;
        else           
            same_energy_count = 0;
            last_total_energy = total_energy;
        end  
    end
    route = zeros(map_length, 2);
    for i = 1 : map_length
        order = find(spin(i, :) == 1, 1);
        route(order, :) = map(i, :);
    end
end
function plot_final_graph(route, width)
    city_num = size(route, 1);
    figure();
    plot(route(:,1),route(:,2),'r.');
    hold on;
    for i = 1 : city_num
        start = route(i,:);
        if i < city_num
            next = route(i+1,:);
        else
            next = route(1,:);
        end
        plot([start(1,1), next(1,1)],[start(1,2), next(1,2)]);
        hold on;
    end
    axis 'equal';
    xlim([width(1, 1), width(2, 1)]);
    ylim([width(1, 2), width(2, 2)]);
    title('Final Graph');
end

