function distance = calculate_distance(start, over, distance_mode)
    switch distance_mode
        case 1
            distance = sqrt(sum((start - over) .^ 2, 2));
        case 2
            distance = sum(abs(start - over), 2);
    end    
end