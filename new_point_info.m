function[]=new_point_info(i,j,d_x,d_y,right_point_x,upper_point_y,nonlocal_size)
    global left_info
    global right_info
    global down_info
    global upper_info
    left_x=nonlocal_size(1);
    right_x=nonlocal_size(2);
    down_y=nonlocal_size(3);
    upper_y=nonlocal_size(4);
    % x direction configuration
    left_point_x=1;
    if((i>left_point_x+left_x)&&(i<right_point_x+1))
        left_info(i,j,3)=0;
        left_info(i,j,2)=left_x;
        left_info(i,j,1)=(i-left_x-1)*d_x;
    elseif((i>left_point_x-1)&&(i<right_point_x+1))
        left_info(i,j,3)=1;
        left_info(i,j,2)=i-left_point_x;
        left_info(i,j,1)=0;
    end
    if((i<right_point_x-right_x)&&(i>left_point_x-1))
        right_info(i,j,3)=0;
        right_info(i,j,2)=right_x;
        right_info(i,j,1)=(i+right_x-1)*d_x;
    elseif((i<right_point_x+1)&&(i>left_point_x-1))
        right_info(i,j,3)=1;
        right_info(i,j,2)=right_point_x-i;
        right_info(i,j,1)=(right_point_x-1)*d_x;
    end 
    %y direction configuration
    down_point_y=1;
    if((j>down_point_y+down_y)&&(j<upper_point_y+1))
        down_info(i,j,3)=0;
        down_info(i,j,2)=down_y;
        down_info(i,j,1)=(j-down_y-1)*d_y;
    elseif((j>down_point_y-1)&&(j<upper_point_y+1))
        down_info(i,j,3)=1;
        down_info(i,j,2)=j-down_point_y;
        down_info(i,j,1)=0;
    end
    if((j<upper_point_y-upper_y)&&(j>down_point_y-1))
        upper_info(i,j,3)=0;
        upper_info(i,j,2)=upper_y;
        upper_info(i,j,1)=(j+upper_y-1)*d_y;
    elseif((j<upper_point_y+1)&&(j>down_point_y-1))
        upper_info(i,j,3)=1;
        upper_info(i,j,2)=upper_point_y-j;
        upper_info(i,j,1)=(upper_point_y-1)*d_y;
    end
end