function[]=new_boun_info(i,global_x,global_y,nonlocal_size)
    global boun_left_info
    global boun_right_info
    global boun_down_info
    global boun_upper_info
    global boun_point_posi
    global right_out_boun_range
    global upper_out_boun_range
    global in_cross_boun_mark
    global d_x
    global d_y
    
    % boun_info(1):end point position
    % boun_info(2):nonlocal range
    in_boun_1=in_cross_boun_mark(global_x,global_y,1);
    in_boun_2=in_cross_boun_mark(global_x,global_y,2);
    left_x=nonlocal_size(1);
    right_x=nonlocal_size(2);
    down_y=nonlocal_size(3);
    upper_y=nonlocal_size(4);
    if(in_boun_1==-1)
        if(in_boun_2==2)
            boun_point_posi(i,1,1)=0;
            boun_point_posi(i,1,2)=global_y*d_y;
            boun_left_info(i,1,1)=0;
            boun_left_info(i,1,2)=0;
            r_range=min(right_out_boun_range(global_y+1),global_x+right_x);
            boun_right_info(i,1,1)=(r_range-1)*d_x;
            boun_right_info(i,1,2)=r_range-global_x;
            d_range=max(global_y+1-down_y,1);
            boun_down_info(i,1,1)=(d_range-1)*d_y;
            boun_down_info(i,1,2)=global_y+1-d_range;
            boun_upper_info(i,1,1)=global_y*d_y;
            boun_upper_info(i,1,2)=0;
        else
            boun_point_posi(i,1,1)=global_x*d_x;
            boun_point_posi(i,1,2)=0;
            l_range=max(global_x+1-left_x,1);
            boun_left_info(i,1,1)=(l_range-1)*d_x;
            boun_left_info(i,1,2)=global_x+1-l_range;
            boun_right_info(i,1,1)=global_x*d_x;
            boun_right_info(i,1,2)=0;
            boun_down_info(i,1,1)=0;
            boun_down_info(i,1,2)=0;
            u_range=min(upper_out_boun_range(global_x+1),global_y+upper_y);
            boun_upper_info(i,1,1)=(u_range-1)*d_y;
            boun_upper_info(i,1,2)=u_range-global_y;
        end
    end
    if((in_boun_1==1)||(in_boun_1==5))
        if(((in_boun_1==1)&&(in_boun_2==1))||((in_boun_1==5)&&(in_boun_2==2)))
            boun_point_posi(i,1,1)=global_x*d_x;
            boun_point_posi(i,1,2)=(global_y-1)*d_y;
            l_range=max(global_x+1-left_x,1);
            boun_left_info(i,1,1)=(l_range-1)*d_x;
            boun_left_info(i,1,2)=global_x+1-l_range;
            r_range=min(global_x+1,global_x+1+right_x);
            boun_right_info(i,1,1)=(r_range-1)*d_x;
            boun_right_info(i,1,2)=r_range-(global_x+1);
            d_range=max(global_y-down_y,1);
            boun_down_info(i,1,1)=(d_range-1)*d_y;
            boun_down_info(i,1,2)=global_y-d_range;
            u_range=min(upper_out_boun_range(global_x+1),global_y+upper_y);
            boun_upper_info(i,1,1)=(u_range-1)*d_y;
            boun_upper_info(i,1,2)=u_range-global_y;
        elseif(((in_boun_1==1)&&(in_boun_2==2))||((in_boun_1==5)&&(in_boun_2==1)))
            boun_point_posi(i,1,1)=(global_x-1)*d_x;
            boun_point_posi(i,1,2)=global_y*d_y;
            l_range=max(global_x-left_x,1);
            boun_left_info(i,1,1)=(l_range-1)*d_x;
            boun_left_info(i,1,2)=global_x-l_range;
            r_range=min(right_out_boun_range(global_y+1),global_x+right_x);
            boun_right_info(i,1,1)=(r_range-1)*d_x;
            boun_right_info(i,1,2)=r_range-global_x;
            d_range=max(global_y+1-down_y,1);
            boun_down_info(i,1,1)=(d_range-1)*d_y;
            boun_down_info(i,1,2)=global_y+1-d_range;
            u_range=min(global_y+1,global_y+1+upper_y);
            boun_upper_info(i,1,1)=(u_range-1)*d_y;
            boun_upper_info(i,1,2)=u_range-(global_y+1);
        end
    elseif((in_boun_1==3)||(in_boun_1==4))
        if(((in_boun_1==3)&&(in_boun_2==1))||((in_boun_1==4)&&((in_boun_2==1)||(in_boun_2==2))))
            boun_point_posi(i,1,1)=global_x*d_x;
            boun_point_posi(i,1,2)=(global_y-1)*d_y;
            l_range=max(global_x+1-left_x,1);
            boun_left_info(i,1,1)=(l_range-1)*d_x;
            boun_left_info(i,1,2)=global_x+1-l_range;
            r_range=min(global_x+1,global_x+1+right_x);
            boun_right_info(i,1,1)=(r_range-1)*d_x;
            boun_right_info(i,1,2)=r_range-(global_x+1);
            d_range=max(global_y-down_y,1);
            boun_down_info(i,1,1)=(d_range-1)*d_y;
            boun_down_info(i,1,2)=global_y-d_range;
            u_range=min(upper_out_boun_range(global_x+1),global_y+upper_y);
            boun_upper_info(i,1,1)=(u_range-1)*d_y;
            boun_upper_info(i,1,2)=u_range-global_y;
            
            boun_point_posi(i,2,1)=(global_x-1)*d_x;
            boun_point_posi(i,2,2)=global_y*d_y;
            l_range=max(global_x-left_x,1);
            boun_left_info(i,2,1)=(l_range-1)*d_x;
            boun_left_info(i,2,2)=global_x-l_range;
            r_range=min(right_out_boun_range(global_y+1),global_x+right_x);
            boun_right_info(i,2,1)=(r_range-1)*d_x;
            boun_right_info(i,2,2)=r_range-global_x;
            d_range=max(global_y+1-down_y,1);
            boun_down_info(i,2,1)=(d_range-1)*d_y;
            boun_down_info(i,2,2)=global_y+1-d_range;
            u_range=min(global_y+1,global_y+1+upper_y);
            boun_upper_info(i,2,1)=(u_range-1)*d_y;
            boun_upper_info(i,2,2)=u_range-(global_y+1);
        end
    end
end