function[]=new_boun_range_info(global_x,global_y)
    global right_out_boun_range
    global upper_out_boun_range
    global point_material
    global out_cross_boun_mark
    
    if(point_material(global_x+1,global_y)==1)
        k=1;
        while((point_material(global_x+k,global_y)==1))
            k=k+1;
        end
        right_out_boun_range(global_y)=global_x+k-1;
        k=1;
        while((point_material(global_x+1,global_y+k)==1))
            k=k+1;
        end
        upper_out_boun_range(global_x+1)=global_y+k-1;
        if(point_material(global_x,global_y+1)==1)
            k=1;
            while((point_material(global_x,global_y+k)==1))
                k=k+1;
            end
            upper_out_boun_range(global_x)=global_y+k-1;
            k=1;
            while((point_material(global_x+k,global_y+1)==1))
                k=k+1;
            end
            right_out_boun_range(global_y+1)=global_x+k-1;
        else
            
        end
    else
        if(point_material(global_x,global_y+1)==1)
            k=1;
            while((point_material(global_x,global_y+k)==1))
                k=k+1;
            end
            upper_out_boun_range(global_x)=global_y+k-1;
            k=1;
            while((point_material(global_x+k,global_y+1)==1))
                k=k+1;
            end
            right_out_boun_range(global_y+1)=global_x+k-1;
        else
            
        end
    end
end