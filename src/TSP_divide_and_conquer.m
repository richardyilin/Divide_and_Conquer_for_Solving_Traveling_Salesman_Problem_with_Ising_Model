close all;
clc;
clear;
% start reading input file
file_name = get_input_file();
fileID = fopen(file_name,'r');
buffer = fscanf(fileID,'%f');
fclose(fileID);
index = 1;
clustering_algorithm = buffer(index, 1); % 1 is kmeans, 2 is agglomerative
index = index + 1;
num_group = buffer(index,1);
index = index + 1;
exponent_start = buffer(index,1);
index = index + 1;    
exponent_increment = buffer(index,1);
index = index + 1;
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
anchor = width(2,:) + 1;
pointer = containers.Map;
[start, final, pointer] =  divide_and_conquer(map, anchor, num_group, distance_mode, clustering_algorithm, pointer, exponent_start, exponent_increment);
final_encoded = encode_key(final);
pointer(final_encoded) = start;
remove(pointer, encode_key(anchor));
plot_final_graph(start, width, pointer, city_num);

function [start, final, pointer] =  divide_and_conquer(map, anchor, num_group, distance_mode, clustering_algorithm, pointer, exponent_start, exponent_increment)
    map_length = size(map, 1);
    if map_length <= num_group
        if map_length > 1
            route = ising_model(map,exponent_start,exponent_increment, distance_mode, anchor);
        else
            route = map;
        end
        start = route(1, :);
        final = route(end, :);
        for i = 1 : map_length - 1
            cur_start = route(i, :);
            cur_start_encoded = encode_key(cur_start);
            cur_final = route(i+1, :);
            pointer(cur_start_encoded) = cur_final;
        end            
    else
        switch clustering_algorithm
            case 1
                switch distance_mode
                    case 1
                        dist_str = 'sqeuclidean';
                    case 2
                        dist_str = 'cityblock';
                end
                [cluster_indices,centroid] = kmeans(map, num_group,'Display','off','Distance',dist_str,...
                    'EmptyAction','error','MaxIter',1000,'OnlinePhase','on','Options',...
                    statset('UseParallel',1),'Replicates',5); % idx is the group each point is assigned to
                centroid_len = length(unique(cluster_indices));
            case 2
                switch distance_mode
                    case 1
                        dist_str = 'euclidean';
                    case 2
                        dist_str = 'cityblock';
                end
                Z = linkage(map,'single',dist_str);
                cluster_indices = cluster(Z,'Maxclust', num_group);
                centroid_len = length(unique(cluster_indices));
                centroid = zeros(centroid_len, 2);
                num_of_points_in_cluster = zeros(centroid_len, 1);
                for i = 1 : map_length
                    cluster_index = cluster_indices(i, 1);
                    centroid(cluster_index, :) = centroid(cluster_index, :) + map(i, :);
                    num_of_points_in_cluster(cluster_index, 1) = num_of_points_in_cluster(cluster_index, 1) + 1;
                end
                centroid = centroid ./ num_of_points_in_cluster; % now is the real centroid
        end
        clusters = cell(centroid_len, 2); % 1 is centroid, 2 is the points in the cluster, cluster{i, 1} is the ith centroid
        centroid_order = ising_model(centroid,exponent_start,exponent_increment, distance_mode, anchor);
        centroid_original_order_to_sorted_order = zeros(centroid_len, 1);
        for i = 1 : centroid_len
            cur_centroid = centroid(i, :);
            sorted_order = find(ismember(centroid_order, cur_centroid, 'row'));
            centroid_original_order_to_sorted_order(i, 1) = sorted_order;
        end
        for i = 1 : centroid_len % add centroid to cell cluster
            sorted_order = centroid_original_order_to_sorted_order(i, 1);
            clusters{sorted_order, 1} = centroid(i, :);
        end
        for i = 1 : map_length % add points to the cluster they belong to
            cluster_index = cluster_indices(i, 1); % the index of the centroid this point belongs to
            sorted_order = centroid_original_order_to_sorted_order(cluster_index, 1);
            clusters{sorted_order, 2}(end+1, :) = map(i, :);
        end
        for i = 1 : centroid_len
            [cluster_start, cluster_final, pointer] =  divide_and_conquer(clusters{i, 2}, anchor, num_group, distance_mode, clustering_algorithm, pointer, exponent_start, exponent_increment);
            if i == 1
                start = cluster_start;
            end
            anchor_encoded = encode_key(anchor);
            pointer(anchor_encoded) = cluster_start;
            anchor = cluster_final;
        end
        final = anchor;
    end
end

function route = ising_model(map,exponent_start,exponent_increment, distance_mode, anchor)
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
    exponent = exponent_start;
    A = mpower(max_weight, exponent);
    init_t = max_weight;
    min_t = (error_rate/(map_length * (log2(map_length) * 3 / 2)));
    beta = 1.0 / ( max_weight * infinite_factor);
    same_energy_count = 0;
    current_t=init_t;
    spin = false(map_length);
    best_spin = false(map_length);
    first_sum_v = map_length; 
    first_sum_j = zeros(map_length,1);
    second_sum_j = map_length;
    second_sum_v = zeros(map_length,1);
    third_sum = 0;
    total_energy = A * (first_sum_v + second_sum_j);
    last_total_energy = Inf;
    best_first_sum_v = map_length; 
    best_first_sum_j = zeros(map_length,1);
    best_second_sum_j = map_length;
    best_second_sum_v = zeros(map_length,1);
    best_third_sum = 0;
    best_total_energy = Inf;
    optimal_spin = false(map_length);
    optimal_total_energy = Inf;
    while (isinf(optimal_total_energy))
        while(current_t>min_t)            
            for v = 1 : map_length
                for j = 1 : map_length
                    if(spin(v,j))% 1 to 0
                        first_sum_j_next = first_sum_j(v,1) - 1;
                        second_sum_v_next = second_sum_v(j,1) - 1; 
                        third_difference = asymmetric_traveling_salesman_calculate_third_difference(spin, data, v, j, false) ;                                           
                    else % 0 to 1
                        first_sum_j_next = first_sum_j(v,1) + 1;
                        second_sum_v_next = second_sum_v(j,1) + 1;
                        third_difference = asymmetric_traveling_salesman_calculate_third_difference(spin, data, v, j, true) ;     
                    end
                    first_sum_v_next = first_sum_v - mpower((1-first_sum_j(v,1)),2) + mpower((1-first_sum_j_next),2);
                    second_sum_j_next = second_sum_j - mpower((1-second_sum_v(j,1)),2) + mpower((1-second_sum_v_next),2);
                    energy_difference = (A * (first_sum_v_next - first_sum_v + second_sum_j_next - second_sum_j)) + third_difference;
                    if (energy_difference < 0)
                        spin(v,j) = (~spin(v,j));
                        first_sum_v = first_sum_v_next; 
                        first_sum_j(v,1) = first_sum_j_next;
                        second_sum_j = second_sum_j_next;
                        second_sum_v(j,1) = second_sum_v_next;
                        third_sum = third_sum + third_difference;
                        total_energy = A * (first_sum_v + second_sum_j) +third_sum;  
                        if (best_total_energy > total_energy)
                            best_spin = spin;
                            best_first_sum_v = first_sum_v; 
                            best_first_sum_j = first_sum_j;
                            best_second_sum_j = second_sum_j;
                            best_second_sum_v = second_sum_v;
                            best_third_sum = third_sum;
                            best_total_energy = total_energy;  
                        end
                        if (optimal_total_energy > total_energy && first_sum_v == 0 && second_sum_j == 0)
                            optimal_spin = spin;
                            optimal_total_energy = total_energy;
                        end
                    else
                        prob = exp(-energy_difference/current_t);
                        if (prob>rand)
                            spin(v,j) = (~spin(v,j));
                            first_sum_v = first_sum_v_next; 
                            first_sum_j(v,1) = first_sum_j_next;
                            second_sum_j = second_sum_j_next;
                            second_sum_v(j,1) = second_sum_v_next;
                            third_sum = third_sum + third_difference;
                            total_energy = A * (first_sum_v + second_sum_j) +third_sum;
                        end                        
                    end      
                end                
            end

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
        if (isinf(optimal_total_energy))
            exponent = exponent + exponent_increment;
            A = mpower(max_weight,exponent);
            same_energy_count = 0;
            current_t=init_t;
            spin = best_spin;
            first_sum_v = best_first_sum_v; 
            first_sum_j = best_first_sum_j;
            second_sum_j = best_second_sum_j;
            second_sum_v = best_second_sum_v;
            third_sum = best_third_sum;
            best_total_energy = A * (best_first_sum_v + best_second_sum_j) + best_third_sum;
            total_energy = best_total_energy;
            last_total_energy = Inf;
        end
    end
    cur_index = find_nearest_point(map, anchor, distance_mode);
    cur_original_order = find(optimal_spin(cur_index,:) == 1,1);
    route = zeros(map_length, 2);
    for real_order = 1 : map_length
        cur_index = find(optimal_spin(:, cur_original_order) == 1,1);
        route(real_order, :) = map(cur_index, :);
        if cur_original_order < map_length
            cur_original_order = cur_original_order + 1;
        else
            cur_original_order = 1;
        end
    end
end

function index_of_nearest_point = find_nearest_point(map, anchor, distance_mode)
    distance = calculate_distance(map, anchor, distance_mode);
    [~, index_of_nearest_point] = min(distance);    
end
function plot_final_graph(start, width, pointer, city_num)
    figure();
    for i = 1 : city_num
        plot(start(1,1),start(1,2),'r.');
        hold on;
        encoded_start = encode_key(start);
        next = pointer(encoded_start);
        plot([start(1,1), next(1,1)],[start(1,2), next(1,2)]);
        hold on;
        start = next;
    end
    axis 'equal';
    xlim([width(1, 1), width(2, 1)]);
    ylim([width(1, 2), width(2, 2)]);
    title('Final Graph');
end

