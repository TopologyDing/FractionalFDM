function[results]=integer_displacement_gradient(i,j,point_type,boun_point_type,max_x,max_y)
    global d_d
    global d_x
    global d_y
    global in_fic_for_out_ux
    global in_fic_for_out_uy
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global in_fic_index
    global out_fic_index
    % compute the integer displacement gradient for each point
    results=sym(zeros(2,1));
    results(1)=0;
    results(2)=0;
    %outer boundary
    if(point_type==1)
        k=in_fic_index(i,j);
        if(boun_point_type==1)
            results(1)=(d_d(i+1,j,1)-in_fic_for_out_ux(k,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-in_fic_for_out_uy(k,1))/(2*d_x);
            results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
            results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
        elseif(boun_point_type==2)
            results(1)=(d_d(i+1,j,1)-d_d(i-1,j,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
            results(3)=(d_d(i,j+1,1)-in_fic_for_out_ux(k,1))/(2*d_y);
            results(4)=(d_d(i,j+1,2)-in_fic_for_out_uy(k,1))/(2*d_y);
        elseif(boun_point_type==3)
            results(1)=(d_d(i+1,j,1)-in_fic_for_out_ux(k,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-in_fic_for_out_uy(k,1))/(2*d_x);
            results(3)=(d_d(i,j+1,1)-in_fic_for_out_ux(k,2))/(2*d_y);
            results(4)=(d_d(i,j+1,2)-in_fic_for_out_uy(k,2))/(2*d_y);
        elseif(boun_point_type==4)
            results(1)=(d_d(i+1,j,1)-d_d(i-1,j,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
            results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
            results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
        else
            if(j==1)
                results(1)=(d_d(i+1,j,1)-in_fic_for_out_ux(k,1))/(2*d_x);
                results(2)=(d_d(i+1,j,2)-in_fic_for_out_uy(k,1))/(2*d_x);
                results(3)=(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y);
                results(4)=(-3*d_d(i,j,2)+4*d_d(i,j+1,2)-d_d(i,j+2,2))/(2*d_y);
            end
            if(i==1)
                results(1)=(-3*d_d(i,j,1)+4*d_d(i+1,j,1)-d_d(i+2,j,1))/(2*d_x);
                results(2)=(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x);
                results(3)=(d_d(i,j+1,1)-in_fic_for_out_ux(k,1))/(2*d_y);
                results(4)=(d_d(i,j+1,2)-in_fic_for_out_uy(k,1))/(2*d_y);
            end
        end
    %inner boundary
    elseif(point_type==-1)
        k=out_fic_index(i,j);
        if(boun_point_type==1)
            if(j~=1)
                results(1)=(out_fic_for_in_ux(k,1)-d_d(i-1,j,1))/(2*d_x);
                results(2)=(out_fic_for_in_uy(k,1)-d_d(i-1,j,2))/(2*d_x);
                results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
            end
        elseif(boun_point_type==2)
            if(i~=1)
                results(1)=(d_d(i+1,j,1)-d_d(i-1,j,1))/(2*d_x);
                results(2)=(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
                results(3)=(out_fic_for_in_ux(k,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(out_fic_for_in_uy(k,1)-d_d(i,j-1,2))/(2*d_y);
            end
        elseif(boun_point_type==3)
            results(1)=(d_d(i+1,j,1)-d_d(i-1,j,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
            results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
            results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
        elseif(boun_point_type==4)
            results(1)=(out_fic_for_in_ux(k,1)-d_d(i-1,j,1))/(2*d_x);
            results(2)=(out_fic_for_in_uy(k,1)-d_d(i-1,j,2))/(2*d_x);
            results(3)=(out_fic_for_in_ux(k,2)-d_d(i,j-1,1))/(2*d_y);
            results(4)=(out_fic_for_in_uy(k,2)-d_d(i,j-1,2))/(2*d_y);
        else
            if(i==1)
                results(1)=(-3*d_d(i,j,1)+4*d_d(i+1,j,1)-d_d(i+2,j,1))/(2*d_x);
                results(2)=(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x);
                results(3)=(out_fic_for_in_ux(k,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(out_fic_for_in_uy(k,1)-d_d(i,j-1,2))/(2*d_y);
            end
            if(j==1)
                results(1)=(out_fic_for_in_ux(k,1)-d_d(i-1,j,1))/(2*d_x);
                results(2)=(out_fic_for_in_uy(k,1)-d_d(i-1,j,2))/(2*d_x);
                results(3)=(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y);
                results(4)=(-3*d_d(i,j,2)+4*d_d(i,j+1,2)-d_d(i,j+2,2))/(2*d_y);
            end
        end
    %points does not involve boundary
    elseif((point_type==2)||(point_type==-2))
        if(i==1)
            results(1)=(-3*d_d(i,j,1)+4*d_d(i+1,j,1)-d_d(i+2,j,1))/(2*d_x);
            results(2)=(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x);
            if(j==1)
                results(3)=(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y);
                results(4)=(-3*d_d(i,j,2)+4*d_d(i,j+1,2)-d_d(i,j+2,2))/(2*d_y);
            elseif(j==max_y)
                results(3)=(3*d_d(i,j,1)-4*d_d(i,j-1,1)+d_d(i,j-2,1))/(2*d_y);
                results(4)=(3*d_d(i,j,2)-4*d_d(i,j-1,2)+d_d(i,j-2,2))/(2*d_y);
            else
                results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
            end
        elseif(i==max_x)
            results(1)=(3*d_d(i,j,1)-4*d_d(i-1,j,1)+d_d(i-2,j,1))/(2*d_x);
            results(2)=(3*d_d(i,j,2)-4*d_d(i-1,j,2)+d_d(i-2,j,2))/(2*d_x);  
            if(j==1)
                results(3)=(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y);
                results(4)=(-3*d_d(i,j,2)+4*d_d(i,j+1,2)-d_d(i,j+2,2))/(2*d_y);
            elseif(j==max_y)
                results(3)=(3*d_d(i,j,1)-4*d_d(i,j-1,1)+d_d(i,j-2,1))/(2*d_y);
                results(4)=(3*d_d(i,j,2)-4*d_d(i,j-1,2)+d_d(i,j-2,2))/(2*d_y);
            else
                results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
            end
        else
            results(1)=(d_d(i+1,j,1)-d_d(i-1,j,1))/(2*d_x);
            results(2)=(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
            if(j==1)
                results(3)=(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y);
                results(4)=(-3*d_d(i,j,2)+4*d_d(i,j+1,2)-d_d(i,j+2,2))/(2*d_y);
            elseif(j==max_y)
                results(3)=(3*d_d(i,j,1)-4*d_d(i,j-1,1)+d_d(i,j-2,1))/(2*d_y);
                results(4)=(3*d_d(i,j,2)-4*d_d(i,j-1,2)+d_d(i,j-2,2))/(2*d_y);
            else
                results(3)=(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
                results(4)=(d_d(i,j+1,2)-d_d(i,j-1,2))/(2*d_y);
            end
        end
    end
    %% Rewrite integer displacement gradient without using symbolic variables
%     switch point_type
%         case 1
%             
%         case -1
%             
%         case 3
%             
%         case 4
%             
%         otherwise
%             
%     end
%     
%     
end