function third_difference = asymmetric_traveling_salesman_calculate_third_difference(spin, data, v, j, plus)
    third_difference = 0;
    len = size(data, 1);
    if(j > 1 && j < len)
        for ver = 1 : len % cancel the punishment j-1 to j           
            if (spin(ver,j+1)) % spin[v][j]  is xuj
                if(plus)
                    third_difference = third_difference + data(v,ver);
                else
                    third_difference = third_difference - data(v,ver);
                end               
            end
            if (spin(ver,j-1)) %spin[v][j]  is x v j+1
                if(plus)
                    third_difference = third_difference + data(ver,v);
                else
                    third_difference = third_difference - data(ver,v);
                end
            end   
        end
    elseif(j == len)
        for ver = 1 : len % cancel the punishment j-1 to j           
            if (spin(ver,1)) % spin[v][j]  is xuj
                if(plus)
                    third_difference = third_difference + data(v,ver);
                else
                    third_difference = third_difference - data(v,ver);
                end               
            end
            if (spin(ver,len-1)) %spin[v][j]  is x v j+1
                if(plus)
                    third_difference = third_difference + data(ver,v);
                else
                    third_difference = third_difference - data(ver,v);
                end
            end   
        end
    else
        for ver = 1 : len % cancel the punishment j-1 to j           
            if (spin(ver,2)) % spin[v][j]  is xuj
                if(plus)
                    third_difference = third_difference + data(v,ver);
                else
                    third_difference = third_difference - data(v,ver);
                end               
            end
            if (spin(ver,len)) %spin[v][j]  is x v j+1
                if(plus)
                    third_difference = third_difference + data(ver,v);
                else
                    third_difference = third_difference - data(ver,v);
                end
            end   
        end
    end
end