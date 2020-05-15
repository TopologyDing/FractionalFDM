function[]=boun_info(in_boun,i,global_x,global_y,nonlocal_size)
    global r
    global d_x
    global d_y
    global boun_left_info
    global boun_right_info
    global boun_down_info
    global boun_upper_info
    global boun_point_posi
    
    in_boun_1=in_boun(1);
    in_boun_2=in_boun(2);
    left_x=nonlocal_size(1);
    down_y=nonlocal_size(3);
    
    left_range=left_x*d_x;
    down_range=down_y*d_y;
    
    a_posi_x=0;
    a_posi_y=0;
    a_posi_dir=0;
    a_left_info_1=0;
    a_left_info_2=0;
    a_left_info_3=0;
    a_right_info_1=0;
    a_right_info_2=0;
    a_right_info_3=0;
    a_down_info_1=0;
    a_down_info_2=0;
    a_down_info_3=0;
    a_upper_info_1=0;
    a_upper_info_2=0;
    a_upper_info_3=0;
    b_posi_dir=0;
    b_posi_x=0;
    b_posi_y=0;
    b_left_info_1=0;
    b_left_info_2=0;
    b_left_info_3=0;
    b_right_info_1=0;
    b_right_info_2=0;
    b_right_info_3=0;
    b_down_info_1=0;
    b_down_info_2=0;
    b_down_info_3=0;
    b_upper_info_1=0;
    b_upper_info_2=0;
    b_upper_info_3=0;
    
    if(in_boun_1==-1)
        if(in_boun_2==2)
            a_posi_x=(global_x-1)*d_x;
            a_posi_dir=2;
            a_posi_y=r;
            a_posi_y_2=0;
            d_range=max(a_posi_y-down_range,a_posi_y_2);
            a_down_info_3=a_posi_y_2>a_posi_y-down_range;
            if(a_posi_y==d_range)
                a_down_info_2=0;
            else
                if(ceil(d_range/d_y)-d_range/d_y<1e-4)
                    b=ceil(d_range/d_y);
                else
                    b=floor(d_range/d_y);
                end
                if((a_posi_y/d_y)-floor(a_posi_y/d_y)<1e-4)
                    a=floor(a_posi_y/d_y);
                else
                    a=ceil(a_posi_y/d_y);
                end
                a_down_info_2=a-b;
            end
            a_down_info_1=d_range;
            a_upper_info_3=0;
            a_upper_info_2=0;
            a_upper_info_1=a_posi_y;
            a_left_info_3=1;
            a_left_info_2=0;
            a_left_info_1=a_posi_x;
            a_right_info_3=0;
            a_right_info_2=0;
            a_right_info_1=a_posi_x;
        else
            a_posi_y=(global_y-1)*d_y;
            a_posi_dir=1;
            a_posi_x=r;
            a_posi_x_2=0;
            l_range=max(a_posi_x-left_range,a_posi_x_2);
            a_left_info_3=a_posi_x_2>a_posi_x-left_range;
            if(a_posi_x==l_range)
                a_left_info_2=0;
            else
                if(ceil(l_range/d_x)-l_range/d_x<1e-4)
                    b=ceil(l_range/d_x);
                else
                    b=floor(l_range/d_x);
                end
                if((a_posi_x/d_x)-floor(a_posi_x/d_x)<1e-4)
                    a=floor(a_posi_x/d_x);
                else
                    a=ceil(a_posi_x/d_x);
                end
                a_left_info_2=a-b;
            end
            a_left_info_1=l_range;
            a_right_info_3=0;
            a_right_info_2=0;
            a_right_info_1=a_posi_x;
            a_upper_info_3=0;
            a_upper_info_2=0;
            a_upper_info_1=a_posi_y;
            a_down_info_3=1;
            a_down_info_2=0;
            a_down_info_1=a_posi_y;
        end
        boun_point_posi(i,1,1)=a_posi_x;
        boun_point_posi(i,1,2)=a_posi_y;
        boun_point_posi(i,1,3)=a_posi_dir;
        boun_left_info(i,1,1)=a_left_info_1;
        boun_left_info(i,1,2)=a_left_info_2;
        boun_left_info(i,1,3)=a_left_info_3;
        boun_right_info(i,1,1)=a_right_info_1;
        boun_right_info(i,1,2)=a_right_info_2;
        boun_right_info(i,1,3)=a_right_info_3;
        boun_down_info(i,1,1)=a_down_info_1;
        boun_down_info(i,1,2)=a_down_info_2;
        boun_down_info(i,1,3)=a_down_info_3;
        boun_upper_info(i,1,1)=a_upper_info_1;
        boun_upper_info(i,1,2)=a_upper_info_2;
        boun_upper_info(i,1,3)=a_upper_info_3;
    end
    if((in_boun_1==1)||(in_boun_1==5))
        if(((in_boun_1==1)&&(in_boun_2==1))||((in_boun_1==5)&&(in_boun_2==2)))
            a_posi_y=(global_y-1)*d_y;
            a_posi_dir=1;
            a_posi_x=sqrt(r^2-(a_posi_y)^2);
            a_posi_x_2=0;
            l_range=max(a_posi_x-left_range,a_posi_x_2);
            a_left_info_3=a_posi_x_2>a_posi_x-left_range;
            if(a_posi_x==l_range)
                a_left_info_2=0;
            else
                if(ceil(l_range/d_x)-l_range/d_x<1e-4)
                    b=ceil(l_range/d_x);
                else
                    b=floor(l_range/d_x);
                end
                if((a_posi_x/d_x)-floor(a_posi_x/d_x)<1e-4)
                    a=floor(a_posi_x/d_x);
                else
                    a=ceil(a_posi_x/d_x);
                end
                a_left_info_2=a-b;
            end
            a_left_info_1=l_range;
            a_right_info_3=1;
            a_right_info_2=0;
            a_right_info_1=a_posi_x;
            a_posi_y_2=0;
            a_upper_info_3=1;
            a_upper_info_2=0;
            a_upper_info_1=a_posi_y;
            down_rg=max(global_y-down_y,1);
            a_down_info_3=1>global_y-down_y;
            a_down_info_2=global_y-down_rg;
            a_down_info_1=(down_rg-1)*d_y;
        elseif(((in_boun_1==1)&&(in_boun_2==2))||((in_boun_1==5)&&(in_boun_2==1)))
            a_posi_x=(global_x-1)*d_x;
            a_posi_dir=2;
            a_posi_y=sqrt(r^2-(a_posi_x)^2);
            a_posi_y_2=0;
            d_range=max(a_posi_y-down_range,a_posi_y_2);
            a_down_info_3=a_posi_y_2>a_posi_y-down_range;
            if(a_posi_y==d_range)
                a_down_info_2=0;
            else
                if(ceil(d_range/d_y)-d_range/d_y<1e-4)
                    b=ceil(d_range/d_y);
                else
                    b=floor(d_range/d_y);
                end
                if((a_posi_y/d_y)-floor(a_posi_y/d_y)<1e-4)
                    a=floor(a_posi_y/d_y);
                else
                    a=ceil(a_posi_y/d_y);
                end
                a_down_info_2=a-b;
            end
            a_down_info_1=d_range;
            a_upper_info_3=1;
            a_upper_info_2=0;
            a_upper_info_1=a_posi_y;
            a_posi_x_2=0;
            a_right_info_3=1;
            a_right_info_2=0;
            a_right_info_1=a_posi_x;
            left_rg=max(global_x-left_x,1);
            a_left_info_3=1>global_x-left_x;
            a_left_info_2=global_x-left_rg;
            a_left_info_1=(left_rg-1)*d_x;
        end
        boun_point_posi(i,1,1)=a_posi_x;
        boun_point_posi(i,1,2)=a_posi_y;
        boun_point_posi(i,1,3)=a_posi_dir;
        boun_left_info(i,1,1)=a_left_info_1;
        boun_left_info(i,1,2)=a_left_info_2;
        boun_left_info(i,1,3)=a_left_info_3;
        boun_right_info(i,1,1)=a_right_info_1;
        boun_right_info(i,1,2)=a_right_info_2;
        boun_right_info(i,1,3)=a_right_info_3;
        boun_down_info(i,1,1)=a_down_info_1;
        boun_down_info(i,1,2)=a_down_info_2;
        boun_down_info(i,1,3)=a_down_info_3;
        boun_upper_info(i,1,1)=a_upper_info_1;
        boun_upper_info(i,1,2)=a_upper_info_2;
        boun_upper_info(i,1,3)=a_upper_info_3;
    elseif((in_boun_1==3)||(in_boun_1==4))
        if(((in_boun_1==3)&&(in_boun_2==1))||((in_boun_1==4)&&((in_boun_2==1)||(in_boun_2==2))))
            %right fic
            a_posi_y=(global_y-1)*d_y;
            a_posi_dir=1;
            a_posi_x=sqrt(r^2-(a_posi_y)^2);
            a_posi_x_2=0;
            l_range=max(a_posi_x-left_range,a_posi_x_2);
            a_left_info_3=a_posi_x_2>a_posi_x-left_range;
            if(a_posi_x==l_range)
                a_left_info_2=0;
            else
                if(ceil(l_range/d_x)-l_range/d_x<1e-4)
                    b=ceil(l_range/d_x);
                else
                    b=floor(l_range/d_x);
                end
                if((a_posi_x/d_x)-floor(a_posi_x/d_x)<1e-4)
                    a=floor(a_posi_x/d_x);
                else
                    a=ceil(a_posi_x/d_x);
                end
                a_left_info_2=a-b;
            end
            a_left_info_1=l_range;
            a_right_info_3=1;
            a_right_info_2=0;
            a_right_info_1=a_posi_x;
            a_posi_y_2=0;
            a_upper_info_3=1;
            a_upper_info_2=0;
            a_upper_info_1=a_posi_y;
            down_rg=max(global_y-down_y,1);
            a_down_info_3=1>global_y-down_y;
            a_down_info_2=global_y-down_rg;
            a_down_info_1=(down_rg-1)*d_y;
            %upper fic
            b_posi_x=(global_x-1)*d_x;
            b_posi_dir=2;
            b_posi_y=sqrt(r^2-(b_posi_x)^2);
            b_posi_y_2=0;
            d_range=max(b_posi_y-down_range,b_posi_y_2);
            b_down_info_3=b_posi_y_2>b_posi_y-down_range;
            if(b_posi_y==d_range)
                b_down_info_2=0;
            else
                if(ceil(d_range/d_y)-d_range/d_y<1e-4)
                    b=ceil(d_range/d_y);
                else
                    b=floor(d_range/d_y);
                end
                if((b_posi_y/d_y)-floor(b_posi_y/d_y)<1e-4)
                    a=floor(b_posi_y/d_y);
                else
                    a=ceil(b_posi_y/d_y);
                end
                b_down_info_2=a-b;
            end
            b_down_info_1=d_range;
            b_upper_info_3=1;
            b_upper_info_2=0;
            b_upper_info_1=b_posi_y;
            b_posi_x_2=0;
            b_right_info_3=1;
            b_right_info_2=0;
            b_right_info_1=b_posi_x;
            left_rg=max(global_x-left_x,1);
            b_left_info_3=1>global_x-left_x;
            b_left_info_2=global_x-left_rg;
            b_left_info_1=(left_rg-1)*d_x;
        end
        boun_point_posi(i,1,1)=a_posi_x;
        boun_point_posi(i,1,2)=a_posi_y;
        boun_point_posi(i,1,3)=a_posi_dir;
        boun_left_info(i,1,1)=a_left_info_1;
        boun_left_info(i,1,2)=a_left_info_2;
        boun_left_info(i,1,3)=a_left_info_3;
        boun_right_info(i,1,1)=a_right_info_1;
        boun_right_info(i,1,2)=a_right_info_2;
        boun_right_info(i,1,3)=a_right_info_3;
        boun_down_info(i,1,1)=a_down_info_1;
        boun_down_info(i,1,2)=a_down_info_2;
        boun_down_info(i,1,3)=a_down_info_3;
        boun_upper_info(i,1,1)=a_upper_info_1;
        boun_upper_info(i,1,2)=a_upper_info_2;
        boun_upper_info(i,1,3)=a_upper_info_3;
        boun_point_posi(i,2,1)=b_posi_x;
        boun_point_posi(i,2,2)=b_posi_y;
        boun_point_posi(i,2,3)=b_posi_dir;
        boun_left_info(i,2,1)=b_left_info_1;
        boun_left_info(i,2,2)=b_left_info_2;
        boun_left_info(i,2,3)=b_left_info_3;
        boun_right_info(i,2,1)=b_right_info_1;
        boun_right_info(i,2,2)=b_right_info_2;
        boun_right_info(i,2,3)=b_right_info_3;
        boun_down_info(i,2,1)=b_down_info_1;
        boun_down_info(i,2,2)=b_down_info_2;
        boun_down_info(i,2,3)=b_down_info_3;
        boun_upper_info(i,2,1)=b_upper_info_1;
        boun_upper_info(i,2,2)=b_upper_info_2;
        boun_upper_info(i,2,3)=b_upper_info_3;
    end
end