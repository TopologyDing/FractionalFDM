function[results]=point_deri_with_one_coord_online(posi,coord,deri_direction,online_dir)
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global out_fic_index
    global out_dir_index
    global point_material
    global d_x
    global d_y
    global d_d
    
    if(online_dir~=deri_direction)
        if(deri_direction==1)
            % deri_direction==1: on y coord
            floor_x=floor(posi/d_x)+1;
            ceil_x=floor_x+1;
            j=coord;
            if(floor_x~=1)
                if(point_material(ceil_x,j)>0)
                    s=out_fic_index(floor_x,j);
                    for k=1:1:2
                        if(out_dir_index(s,k)==1)
                            right_x=out_fic_for_in_ux(s,k);
                            right_y=out_fic_for_in_uy(s,k);
                            break;
                        end
                    end
                    x_disp=[d_d(floor_x-1,j,1);d_d(floor_x,j,1);right_x];
                    y_disp=[d_d(floor_x-1,j,2);d_d(floor_x,j,2);right_y];
                    results(1)=second_order_interpolation(x_disp,posi,floor_x,1,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,floor_x,1,d_x,d_y,2);
                else
                    x_disp=[d_d(floor_x-1,j,1);d_d(floor_x,j,1);d_d(ceil_x,j,1)];
                    y_disp=[d_d(floor_x-1,j,2);d_d(floor_x,j,2);d_d(ceil_x,j,2)];
                    results(1)=second_order_interpolation(x_disp,posi,floor_x,1,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,floor_x,1,d_x,d_y,2);
                end
            else
                if(point_material(ceil_x+1,j)>0)
                    s=out_fic_index(ceil_x,j);
                    for k=1:1:2
                        if(out_dir_index(s,k)==1)
                            right_x=out_fic_for_in_ux(s,k);
                            right_y=out_fic_for_in_uy(s,k);
                            break;
                        end
                    end
                    x_disp=[d_d(floor_x,j,1);d_d(ceil_x,j,1);right_x];
                    y_disp=[d_d(floor_x,j,2);d_d(ceil_x,j,2);right_y];
                    results(1)=second_order_interpolation(x_disp,posi,ceil_x,1,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,ceil_x,1,d_x,d_y,2);
                else
                    x_disp=[d_d(floor_x,j,1);d_d(ceil_x,j,1);d_d(ceil_x+1,j,1)];
                    y_disp=[d_d(floor_x,j,2);d_d(ceil_x,j,2);d_d(ceil_x+1,j,2)];
                    results(1)=second_order_interpolation(x_disp,posi,ceil_x,1,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,ceil_x,1,d_x,d_y,2);
                end
            end
        elseif(deri_direction==2)
            % deri_direction==2: on x coord
            floor_y=floor(posi/d_y)+1;
            ceil_y=floor_y+1;
            i=coord;
            if(floor_y~=1)
                if(point_material(i,ceil_y)>0)
                    s=out_fic_index(i,floor_y);
                    for k=1:1:2
                        if(out_dir_index(s,k)==2)
                            upper_x=out_fic_for_in_ux(s,k);
                            upper_y=out_fic_for_in_uy(s,k);
                            break;
                        end
                    end
                    x_disp=[d_d(i,floor_y-1,1);d_d(i,floor_y,1);upper_x];
                    y_disp=[d_d(i,floor_y-1,2);d_d(i,floor_y,2);upper_y];
                    results(1)=second_order_interpolation(x_disp,posi,floor_y,0,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,floor_y,0,d_x,d_y,2);
                else
                    x_disp=[d_d(i,floor_y-1,1);d_d(i,floor_y,1);d_d(i,ceil_y,1)];
                    y_disp=[d_d(i,floor_y-1,2);d_d(i,floor_y,2);d_d(i,ceil_y,2)];
                    results(1)=second_order_interpolation(x_disp,posi,floor_y,0,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,floor_y,0,d_x,d_y,2);
                end
            else
                if(point_material(i,ceil_y+1)>0)
                    s=out_fic_index(i,ceil_y);
                    for k=1:1:2
                        if(out_dir_index(s,k)==2)
                            upper_x=out_fic_for_in_ux(s,k);
                            upper_y=out_fic_for_in_uy(s,k);
                            break;
                        end
                    end
                    x_disp=[d_d(i,floor_y,1);d_d(i,ceil_y,1);upper_x];
                    y_disp=[d_d(i,floor_y,2);d_d(i,ceil_y,2);upper_y];
                    results(1)=second_order_interpolation(x_disp,posi,ceil_y,0,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,ceil_y,0,d_x,d_y,2);
                else
                    x_disp=[d_d(i,floor_y,1);d_d(i,ceil_y,1);d_d(i,ceil_y+1,1)];
                    y_disp=[d_d(i,floor_y,2);d_d(i,ceil_y,2);d_d(i,ceil_y+1,2)];
                    results(1)=second_order_interpolation(x_disp,posi,ceil_y,0,d_x,d_y,2);
                    results(2)=second_order_interpolation(y_disp,posi,ceil_y,0,d_x,d_y,2);
                end
            end
        end
    else
        % this condition can be computed by random point interpolation
    end
end