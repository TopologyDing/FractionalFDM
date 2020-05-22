function[results]=random_point_derivative(posi,deri_direction)
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global out_fic_index
    global out_dir_index
    global in_cross_boun_mark
    global point_material
    global d_x
    global d_y
    global d_d
    % deri_direction==1:
    % compute x direction derivative
    % deri_direction==2:
    % compute y direction derivative
    posi_x=posi(1);
    posi_y=posi(2);
    floor_x=floor(posi_x/d_x)+1;
    ceil_x=floor_x+1;
    floor_y=floor(posi_y/d_y)+1;
    ceil_y=floor_y+1;
    
    top_right_point=point_material(ceil_x+1-(floor_x>1),ceil_y+1-(floor_y>1));
    if(top_right_point<0)
        if((floor_x>1)&&(floor_y>1))
            center_coord(1)=floor_x;
            center_coord(2)=floor_y;
            nine_point_x=[d_d(floor_x-1,ceil_y,1) d_d(floor_x,ceil_y,1) d_d(ceil_x,ceil_y,1);
            d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1);
            d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
            nine_point_y=[d_d(floor_x-1,ceil_y,2) d_d(floor_x,ceil_y,2) d_d(ceil_x,ceil_y,2);
            d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2);
            d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
            results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
        elseif((floor_x==1)&&(floor_y>1))
            center_coord(1)=ceil_x;
            center_coord(2)=floor_y;
            nine_point_x=[d_d(floor_x,ceil_y,1) d_d(floor_x+1,ceil_y,1) d_d(ceil_x+1,ceil_y,1);
            d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1);
            d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
            nine_point_y=[d_d(floor_x,ceil_y,2) d_d(floor_x+1,ceil_y,2) d_d(ceil_x+1,ceil_y,2);
            d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2);
            d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
            results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
        elseif((floor_x==1)&&(floor_y==1))
            center_coord(1)=ceil_x;
            center_coord(2)=ceil_y;
            nine_point_x=[d_d(floor_x,ceil_y+1,1) d_d(floor_x+1,ceil_y+1,1) d_d(ceil_x+1,ceil_y+1,1);
            d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_d(ceil_x+1,floor_y+1,1);
            d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
            nine_point_y=[d_d(floor_x,ceil_y+1,2) d_d(floor_x+1,ceil_y+1,2) d_d(ceil_x+1,ceil_y+1,2);
            d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_d(ceil_x+1,floor_y+1,2);
            d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
            results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
        elseif((floor_x>1)&&(floor_y==1))
            center_coord(1)=floor_x;
            center_coord(2)=ceil_y;
            nine_point_x=[d_d(floor_x-1,ceil_y+1,1) d_d(floor_x,ceil_y+1,1) d_d(ceil_x,ceil_y+1,1);
            d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_d(ceil_x,floor_y+1,1);
            d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
            nine_point_y=[d_d(floor_x-1,ceil_y+1,2) d_d(floor_x,ceil_y+1,2) d_d(ceil_x,ceil_y+1,2);
            d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_d(ceil_x,floor_y+1,2);
            d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
            results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
        end 
    else
        if((floor_x>1)&&(floor_y>1))
            center_coord(1)=floor_x;
            center_coord(2)=floor_y;
            if(in_cross_boun_mark(floor_x,floor_y,1)==1)
                if(in_cross_boun_mark(floor_x,floor_y,2)==1)
                    p_1=out_fic_index(floor_x,ceil_y);
                    p_2=out_fic_index(floor_x,floor_y);
                    p_3=out_fic_index(floor_x,floor_y-1);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==1)
                            d2_r_fic_x=out_fic_for_in_ux(p_3,k);
                            d2_r_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p_2,1);
                    d_r_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y,1) d_d(floor_x,ceil_y,1) u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x-1,ceil_y,2) d_d(floor_x,ceil_y,2) u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                elseif(in_cross_boun_mark(floor_x,floor_y,2)==2)
                    p_1=out_fic_index(ceil_x,floor_y);
                    p_2=out_fic_index(floor_x,floor_y);
                    p_3=out_fic_index(floor_x-1,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==2)
                            u_l2_fic_x=out_fic_for_in_ux(p_3,k);
                            u_l2_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p_2,1);
                    u_l_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1);
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2);
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                end
            elseif(in_cross_boun_mark(floor_x,floor_y,1)==2)
                % might consider to use another way
                if(deri_direction==1)
                    p=out_fic_index(ceil_x,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                else
                    p=out_fic_index(floor_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                end
                nine_point_x=[d_d(floor_x-1,ceil_y,1) d_d(floor_x,ceil_y,1) u_r_fic_x;
                d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1);
                d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                nine_point_y=[d_d(floor_x-1,ceil_y,2) d_d(floor_x,ceil_y,2) u_r_fic_y;
                d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2);
                d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
            elseif(in_cross_boun_mark(floor_x,floor_y,1)==3)
                p=out_fic_index(floor_x,floor_y);
                u_r_fic_x=out_fic_for_in_ux(p,3);
                u_r_fic_y=out_fic_for_in_uy(p,3);
                u_l_fic_x=out_fic_for_in_ux(p,2);
                u_l_fic_y=out_fic_for_in_uy(p,2);
                d_r_fic_x=out_fic_for_in_ux(p,1);
                d_r_fic_y=out_fic_for_in_uy(p,1);
                nine_point_x=[d_d(floor_x-1,ceil_y,1) u_l_fic_x u_r_fic_x;
                d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_r_fic_x;
                d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                nine_point_y=[d_d(floor_x-1,ceil_y,2) u_l_fic_y u_r_fic_y;
                d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_r_fic_y;
                d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
            elseif(in_cross_boun_mark(floor_x,floor_y,1)==4)
                p=out_fic_index(floor_x,floor_y);
                if(in_cross_boun_mark(floor_x,floor_y,2)==1)
                    p_1=out_fic_index(floor_x-1,floor_y);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    u_l2_fic_x=out_fic_for_in_ux(p_1,1);
                    u_l2_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                elseif(in_cross_boun_mark(floor_x,floor_y,2)==2)
                    p_1=out_fic_index(floor_x,floor_y-1);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    d2_r_fic_x=out_fic_for_in_ux(p_1,1);
                    d2_r_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x-1,ceil_y,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                end
            elseif(in_cross_boun_mark(floor_x,floor_y,1)==5)
                p=out_fic_index(floor_x,floor_y);
                if(in_cross_boun_mark(floor_x,floor_y,2)==1)
                    p_1=out_fic_index(ceil_x,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p,1);
                    u_l_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1);
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                    nine_point_y=[d_d(floor_x-1,ceil_y,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2);
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                elseif(in_cross_boun_mark(floor_x,floor_y,2)==2)
                    p_1=out_fic_index(floor_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y,1) d_d(floor_x,ceil_y,1) u_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y-1,1) d_d(floor_x,floor_y-1,1) d_d(ceil_x,floor_y-1,1)];
                    nine_point_y=[d_d(floor_x-1,ceil_y,2) d_d(floor_x,ceil_y,2) u_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y-1,2) d_d(floor_x,floor_y-1,2) d_d(ceil_x,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,1,deri_direction);
                end
            end
        elseif((floor_x==1)&&(floor_y>1))
            center_coord(1)=ceil_x;
            center_coord(2)=floor_y;
            if(in_cross_boun_mark(ceil_x,floor_y,1)==1)
                if(in_cross_boun_mark(ceil_x,floor_y,2)==1)
                    p_1=out_fic_index(ceil_x,ceil_y);
                    p_2=out_fic_index(ceil_x,floor_y);
                    p_3=out_fic_index(ceil_x,floor_y-1);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==1)
                            d2_r_fic_x=out_fic_for_in_ux(p_3,k);
                            d2_r_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p_2,1);
                    d_r_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[d_d(ceil_x-1,ceil_y,1) d_d(ceil_x,ceil_y,1) u_r_fic_x;
                    d_d(ceil_x-1,floor_y,1) d_d(ceil_x,floor_y,1) d_r_fic_x;
                    d_d(ceil_x-1,floor_y-1,1) d_d(ceil_x,floor_y-1,1) d2_r_fic_x];
                    nine_point_y=[d_d(ceil_x-1,ceil_y,2) d_d(ceil_x,ceil_y,2) u_r_fic_y;
                    d_d(ceil_x-1,floor_y,2) d_d(ceil_x,floor_y,2) d_r_fic_y;
                    d_d(ceil_x-1,floor_y-1,2) d_d(ceil_x,floor_y-1,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,floor_y,2)==2)
                    p_1=out_fic_index(ceil_x+1,floor_y);
                    p_2=out_fic_index(ceil_x,floor_y);
                    p_3=out_fic_index(floor_x,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==2)
                            u_l2_fic_x=out_fic_for_in_ux(p_3,k);
                            u_l2_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p_2,1);
                    u_l_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1);
                    d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2);
                    d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                end
            elseif(in_cross_boun_mark(ceil_x,floor_y,1)==2)
                if(deri_direction==1)
                    p=out_fic_index(ceil_x+1,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                else
                    p=out_fic_index(ceil_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                end
                nine_point_x=[d_d(floor_x,ceil_y,1) d_d(floor_x+1,ceil_y,1) u_r_fic_x;
                d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1);
                d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                nine_point_y=[d_d(floor_x,ceil_y,2) d_d(floor_x+1,ceil_y,2) u_r_fic_y;
                d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2);
                d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
            elseif(in_cross_boun_mark(ceil_x,floor_y,1)==3)
                p=out_fic_index(ceil_x,floor_y);
                u_r_fic_x=out_fic_for_in_ux(p,3);
                u_r_fic_y=out_fic_for_in_uy(p,3);
                u_l_fic_x=out_fic_for_in_ux(p,2);
                u_l_fic_y=out_fic_for_in_uy(p,2);
                d_r_fic_x=out_fic_for_in_ux(p,1);
                d_r_fic_y=out_fic_for_in_uy(p,1);
                nine_point_x=[d_d(floor_x,ceil_y,1) u_l_fic_x u_r_fic_x;
                d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_r_fic_x;
                d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                nine_point_y=[d_d(floor_x,ceil_y,2) u_l_fic_y u_r_fic_y;
                d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_r_fic_y;
                d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
            elseif(in_cross_boun_mark(ceil_x,floor_y,1)==4)
                p=out_fic_index(ceil_x,floor_y);
                if(in_cross_boun_mark(ceil_x,floor_y,2)==1)
                    p_1=out_fic_index(floor_x,floor_y);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    u_l2_fic_x=out_fic_for_in_ux(p_1,1);
                    u_l2_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_r_fic_x;
                    d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_r_fic_y;
                    d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,floor_y,2)==2)
                    p_1=out_fic_index(ceil_x,floor_y-1);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    d2_r_fic_x=out_fic_for_in_ux(p_1,1);
                    d2_r_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[d_d(floor_x,ceil_y,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_r_fic_x;
                    d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x,ceil_y,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_r_fic_y;
                    d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                end
            elseif(in_cross_boun_mark(ceil_x,floor_y,1)==5)
                p=out_fic_index(ceil_x,floor_y);
                if(in_cross_boun_mark(ceil_x,floor_y,2)==1)
                    p_1=out_fic_index(ceil_x+1,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        break;
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p,1);
                    u_l_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x,ceil_y,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1);
                    d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                    nine_point_y=[d_d(floor_x,ceil_y,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2);
                    d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,floor_y,2)==2)
                    p_1=out_fic_index(ceil_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=otu_fic_for_in_ux(p,1);
                    nine_point_x=[d_d(floor_x,ceil_y,1) d_d(floor_x+1,ceil_y,1) u_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_r_fic_x;
                    d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y-1,1) d_d(ceil_x+1,floor_y-1,1)];
                    nine_point_y=[d_d(floor_x,ceil_y,2) d_d(floor_x+1,ceil_y,2) u_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_r_fic_y;
                    d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y-1,2) d_d(ceil_x+1,floor_y-1,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,2,deri_direction);
                end
            end
        elseif((floor_x==1)&&(floor_y==1))
            center_coord(1)=ceil_x;
            center_coord(2)=ceil_y;
            if(in_cross_boun_mark(ceil_x,ceil_y,1)==1)
                if(in_cross_boun_mark(ceil_x,ceil_y,2)==1)
                    p_1=out_fic_index(ceil_x,ceil_y+1);
                    p_2=out_fic_index(ceil_x,ceil_y);
                    p_3=out_fic_index(ceil_x,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==1)
                            d2_r_fic_x=out_fic_for_in_ux(p_3,k);
                            d2_r_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p_2,1);
                    d_r_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[d_d(ceil_x-1,ceil_y+1,1) d_d(ceil_x,ceil_y+1,1) u_r_fic_x;
                    d_d(ceil_x-1,floor_y+1,1) d_d(ceil_x,floor_y+1,1) d_r_fic_x;
                    d_d(ceil_x-1,floor_y,1) d_d(ceil_x,floor_y,1) d2_r_fic_x];
                    nine_point_y=[d_d(ceil_x-1,ceil_y+1,2) d_d(ceil_x,ceil_y+1,2) u_r_fic_y;
                    d_d(ceil_x-1,floor_y+1,2) d_d(ceil_x,floor_y+1,2) d_r_fic_y;
                    d_d(ceil_x-1,floor_y,2) d_d(ceil_x,floor_y,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,ceil_y,2)==2)
                    p_1=out_fic_index(ceil_x+1,ceil_y);
                    p_2=out_fic_index(ceil_x,ceil_y);
                    p_3=out_fic_index(floor_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==2)
                            u_l2_fic_x=out_fic_for_in_ux(p_3,k);
                            u_l2_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p_2,1);
                    u_l_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_d(ceil_x+1,floor_y+1,1);
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_d(ceil_x+1,floor_y+1,2);
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                end
            elseif(in_cross_boun_mark(ceil_x,ceil_y,1)==2)
                if(deri_direction==1)
                    p=out_fic_index(ceil_x+1,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                else
                    p=out_fic_index(ceil_x,ceil_y+1);
                    for k=1:1:2
                        if(out_dir_index(p,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                end
                nine_point_x=[d_d(floor_x,ceil_y+1,1) d_d(floor_x+1,ceil_y+1,1) u_r_fic_x;
                d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y+1,1) d_d(ceil_x+1,floor_y+1,1);
                d_d(floor_x,floor_y-1,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                nine_point_y=[d_d(floor_x,ceil_y+1,2) d_d(floor_x+1,ceil_y+1,2) u_r_fic_y;
                d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y+1,2) d_d(ceil_x+1,floor_y+1,2);
                d_d(floor_x,floor_y-1,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
            elseif(in_cross_boun_mark(ceil_x,ceil_y,1)==3)
                p=out_fci_index(ceil_x,ceil_y);
                u_r_fic_x=out_fic_for_in_ux(p,3);
                u_r_fic_y=out_fic_for_in_uy(p,3);
                u_l_fic_x=out_fic_for_in_ux(p,2);
                u_l_fic_y=out_fic_for_in_uy(p,2);
                d_r_fic_x=out_fic_for_in_ux(p,1);
                d_r_fic_y=out_fic_for_in_uy(p,1);
                nine_point_x=[d_d(floor_x,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_r_fic_x;
                d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                nine_point_y=[d_d(floor_x,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_r_fic_y;
                d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
            elseif(in_cross_boun_mark(ceil_x,ceil_y,1)==4)
                p=out_fic_index(ceil_x,ceil_y);
                if(in_cross_boun_mark(ceil_x,ceil_y,2)==1)
                    p_1=out_fic_index(floor_x,ceil_y);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    u_l2_fic_x=out_fic_for_in_ux(p_1,1);
                    u_l2_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,ceil_y,2)==2)
                    p_1=out_fic_index(ceil_x+1,floor_y);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    d2_r_fic_x=out_fic_for_in_ux(p_1,1);
                    d2_r_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[d_d(floor_x,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                end
            elseif(in_cross_boun_mark(ceil_x,ceil_y,1)==5)
                p=out_fic_index(ceil_x,ceil_y);
                if(in_cross_boun_mark(ceil_x,ceil_y,2)==1)
                    p_1=out_fic_index(ceil_x+1,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p,1);
                    u_l_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_d(ceil_x+1,floor_y+1,1);
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                    nine_point_y=[d_d(floor_x,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_d(ceil_x+1,floor_y+1,2);
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                elseif(in_cross_boun_mark(ceil_x,ceil_y,2)==2)
                    p_1=out_fic_index(ceil_x,ceil_y+1);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=otu_fic_for_in_ux(p,1);
                    nine_point_x=[d_d(floor_x,ceil_y+1,1) d_d(floor_x+1,ceil_y+1,1) u_r_fic_x;
                    d_d(floor_x,floor_y+1,1) d_d(floor_x+1,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x,floor_y,1) d_d(floor_x+1,floor_y,1) d_d(ceil_x+1,floor_y,1)];
                    nine_point_y=[d_d(floor_x,ceil_y+1,2) d_d(floor_x+1,ceil_y+1,2) u_r_fic_y;
                    d_d(floor_x,floor_y+1,2) d_d(floor_x+1,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x,floor_y,2) d_d(floor_x+1,floor_y,2) d_d(ceil_x+1,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,3,deri_direction);
                end
            end
        elseif((floor_x>1)&&(floor_y==1))
            center_coord(1)=floor_x;
            center_coord(2)=ceil_y;
            if(in_cross_boun_mark(floor_x,ceil_y,1)==1)
                if(in_cross_boun_mark(floor_x,ceil_y,2)==1)
                    p_1=out_fic_index(floor_x,ceil_y+1);
                    p_2=out_fic_index(floor_x,ceil_y);
                    p_3=out_fic_index(floor_x,floor_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==1)
                            d2_r_fic_x=out_fic_for_in_ux(p_3,k);
                            d2_r_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p_2,1);
                    d_r_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y+1,1) d_d(floor_x,ceil_y+1,1) u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x-1,ceil_y+1,2) d_d(floor_x,ceil_y+1,2) u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                elseif(in_cross_boun_mark(floor_x,ceil_y,2)==2)
                    p_1=out_fic_index(ceil_x,ceil_y);
                    p_2=out_fic_index(floor_x,ceil_y);
                    p_3=out_fic_index(floor_x-1,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                        end
                        if(out_dir_index(p_3,k)==2)
                            u_l2_fic_x=out_fic_for_in_ux(p_3,k);
                            u_l2_fic_y=out_fic_for_in_uy(p_3,k);
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p_2,1);
                    u_l_fic_y=out_fic_for_in_uy(p_2,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_d(ceil_x,floor_y+1,1);
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_d(ceil_x,floor_y+1,2);
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                end
            elseif(in_cross_boun_mark(floor_x,ceil_y,1)==2)
                if(deri_direction==1)
                    p=out_fic_index(ceil_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                else
                    p=out_fic_index(floor_x,ceil_y+1);
                    for k=1:1:2
                        if(out_dir_index(p,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p,k);
                            u_r_fic_y=out_fic_for_in_uy(p,k); 
                            break;
                        end
                    end
                end
                nine_point_x=[d_d(floor_x-1,ceil_y+1,1) d_d(floor_x,ceil_y+1,1) u_r_fic_x;
                d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_d(ceil_x,floor_y+1,1);
                d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                nine_point_y=[d_d(floor_x-1,ceil_y+1,2) d_d(floor_x,ceil_y+1,2) u_r_fic_y;
                d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_d(ceil_x,floor_y+1,2);
                d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
            elseif(in_cross_boun_mark(floor_x,ceil_y,1)==3)
                p=out_fic_index(floor_x,ceil_y);
                u_r_fic_x=out_fic_for_in_ux(p,3);
                u_r_fic_y=out_fic_for_in_uy(p,3);
                u_l_fic_x=out_fic_for_in_ux(p,2);
                u_l_fic_y=out_fic_for_in_uy(p,2);
                d_r_fic_x=out_fic_for_in_ux(p,1);
                d_r_fic_y=out_fic_for_in_uy(p,1);
                nine_point_x=[d_d(floor_x-1,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_r_fic_x;
                d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                nine_point_y=[d_d(floor_x-1,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_r_fic_y;
                d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
            elseif(in_cross_boun_mark(floor_x,ceil_y,1)==4)
                p=out_fic_index(floor_x,ceil_y);
                if(in_cross_boun_mark(floor_x,ceil_y,2)==1)
                    p_1=out_fic_index(floor_x-1,ceil_y);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    u_l2_fic_x=out_fic_for_in_ux(p_1,1);
                    u_l2_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[u_l2_fic_x u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                    nine_point_y=[u_l2_fic_y u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                elseif(in_cross_boun_mark(floor_x,ceil_y,2)==2)
                    p_1=out_fic_index(floor_x,ceil_y-1);
                    u_r_fic_x=out_fic_for_in_ux(p,3);
                    u_r_fic_y=out_fic_for_in_uy(p,3);
                    u_l_fic_x=out_fic_for_in_ux(p,2);
                    u_l_fic_y=out_fic_for_in_uy(p,2);
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    d2_r_fic_x=out_fic_for_in_ux(p_1,1);
                    d2_r_fic_y=out_fic_for_in_uy(p_1,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d2_r_fic_x];
                    nine_point_y=[d_d(floor_x-1,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d2_r_fic_y];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                end
            elseif(in_cross_boun_mark(floor_x,ceil_y,1)==5)
                p=out_fic_index(floor_x,ceil_y);
                if(in_cross_boun_mark(floor_x,ceil_y,2)==1)
                    p_1=out_fic_index(ceil_x,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==2)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    u_l_fic_x=out_fic_for_in_ux(p,1);
                    u_l_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y+1,1) u_l_fic_x u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_d(ceil_x,floor_y+1,1);
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                    nine_point_y=[d_d(floor_x-1,ceil_y+1,2) u_l_fic_y u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_d(ceil_x,floor_y+1,2);
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                elseif(in_cross_boun_mark(floor_x,ceil_y,2)==2)
                    p_1=out_fic_index(floor_x,ceil_y+1);
                    for k=1:1:2
                        if(out_dir_index(p_1,k)==1)
                            u_r_fic_x=out_fic_for_in_ux(p_1,k);
                            u_r_fic_y=out_fic_for_in_uy(p_1,k);
                            break;
                        end
                    end
                    d_r_fic_x=out_fic_for_in_ux(p,1);
                    d_r_fic_y=out_fic_for_in_uy(p,1);
                    nine_point_x=[d_d(floor_x-1,ceil_y+1,1) d_d(floor_x,ceil_y+1,1) u_r_fic_x;
                    d_d(floor_x-1,floor_y+1,1) d_d(floor_x,floor_y+1,1) d_r_fic_x;
                    d_d(floor_x-1,floor_y,1) d_d(floor_x,floor_y,1) d_d(ceil_x,floor_y,1)];
                    nine_point_y=[d_d(floor_x-1,ceil_y+1,2) d_d(floor_x,ceil_y+1,2) u_r_fic_y;
                    d_d(floor_x-1,floor_y+1,2) d_d(floor_x,floor_y+1,2) d_r_fic_y;
                    d_d(floor_x-1,floor_y,2) d_d(floor_x,floor_y,2) d_d(ceil_x,floor_y,2)];
                    results=single_point_derivative(nine_point_x,nine_point_y,posi,center_coord,d_x,d_y,4,deri_direction);
                end
            end    
        end
    end
        
        

end