function[]=boun_derivative_with_nine_points(i,j,k,fic_point_num)
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global in_fic_for_out_ux
    global in_for_for_out_uy
    global in_fic_index
    global in_dir_index
    global out_fic_index
    global out_dir_index
    global d_x
    global d_y
    global d_d
    
    results=sym(zeros(4,1));
    if(in_fic_for_out_ux(k,2)==0)
        if(out_fic_index(i,j,1)==1)
            
        elseif(out_fic_index(i,j,1)==2)
            
        end
    else
        
    end
end