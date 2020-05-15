function[]=overall_boun_equilibrium(high_x,high_y,max_x,max_y,miu_aluminum,miu_steel,lambda_steel,d_x,d_y,size_in_fic_for_out,size_out_fic_for_in)
    global equilibrium
    global sec_order_frac_deri_x_disp_y
    global sec_order_frac_deri_y_disp_x
    global in_fic_for_out_uy
    global in_fic_for_out_ux
    global d_d
    

    for i=1
        for j=1
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=0;
            d_d(i,j,1)=0;
            d_d(i,j,2)=0;
        end
    end
    
    for i=1
        for j=2:1:high_y-1
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=miu_aluminum*(sec_order_frac_deri_x_disp_y(i,j)+sec_order_frac_deri_y_disp_x(i,j));
            d_d(i,j,1)=0;
        end
    end
    %}
    
    for i=1
        for j=high_y
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=miu_steel*((-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                (d_d(i,j+1,1)-in_fic_for_out_ux(1,1))/(2*d_y));
            d_d(i,j,1)=0;
        end
    end
    
    for i=1
        for j=high_y+1:1:max_y-1
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=miu_steel*(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                miu_steel*(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
            d_d(i,j,1)=0;
        end
    end

    for i=2:1:high_x-1
        for j=1
            equilibrium(i,j,1)=miu_aluminum*(sec_order_frac_deri_x_disp_y(i,j)+sec_order_frac_deri_y_disp_x(i,j));
            equilibrium(i,j,2)=0;
            d_d(i,j,2)=0;
        end
    end
    %}
    
    for i=high_x
        for j=1
            equilibrium(i,j,1)=miu_steel*((d_d(i+1,j,2)-in_fic_for_out_uy(size_in_fic_for_out,1))/(2*d_x)+...
                (-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y));
            equilibrium(i,j,2)=0;
            d_d(i,j,2)=0;
        end
    end
            
    for i=high_x+1:1:max_x-1
        for j=1
            equilibrium(i,j,1)=miu_steel*((d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x)+...
                (-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y));
            equilibrium(i,j,2)=0;
            d_d(i,j,2)=0;
        end
    end
    %{
    for i=max_x
        for j=1
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
                lambda_steel*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
            equilibrium(i,j,2)=0;
            d_d(i,j,2)=0;
        end
    end
    
    for i=max_x
        for j=2:1:max_y-1
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,j,1)-4*d_d(max_x-1,j,1)+d_d(max_x-2,j,1))/(2*d_x)+...
                lambda_steel*(d_d(max_x,j+1,2)-d_d(max_x,j-1,2))/(2*d_y);
            equilibrium(i,j,2)=miu_steel*(d_d(max_x,j+1,1)-d_d(max_x,j-1,1))/(2*d_y)+...
                miu_steel*(3*d_d(max_x,j,2)-4*d_d(max_x-1,j,2)+d_d(max_x-2,j,2))/(2*d_x);
        end
    end
    
    for i=1
        for j=max_y
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=lambda_steel*(-3*d_d(1,max_y,1)+4*d_d(2,max_y,1)-d_d(3,max_y,1))/(2*d_x)+...
                (2*miu_steel+lambda_steel)*(3*d_d(1,max_y,2)-4*d_d(1,max_y-1,2)+d_d(1,max_y-2,2))/(2*d_y);
            d_d(i,j,1)=0;
        end
    end
    
    for i=2:1:max_x-1
        for j=max_y
            equilibrium(i,j,2)=lambda_steel*(d_d(i+1,max_y,1)-d_d(i-1,max_y,1))/(2*d_x)+...
                (2*miu_steel+lambda_steel)*(3*d_d(i,max_y,2)-4*d_d(i,max_y-1,2)+d_d(i,max_y-2,2))/(2*d_y);
            equilibrium(i,j,1)=miu_steel*(3*d_d(i,max_y,1)-4*d_d(i,max_y-1,1)+d_d(i,max_y-2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,max_y,2)-d_d(i-1,max_y,2))/(2*d_x);
        end
    end
    
    for i=max_x
        for j=max_y
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
                lambda_steel*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
            equilibrium(i,j,2)=lambda_steel*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
                (2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
        end
    end
    %}
    
    for i=1
        for j=max_y
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=lambda_steel*(-3*d_d(i,j,1)+4*d_d(i+1,j,1)-d_d(i+2,j,1))/(2*d_x)+...
                    (2*miu_steel+lambda_steel)*(3*d_d(i,j,2)-4*d_d(i,j-1,2)+d_d(i,j-2,2))/(2*d_y);
            d_d(i,j,1)=0;
        end
    end
    
    for i=max_x
        for j=1
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
                lambda_steel*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
            equilibrium(i,j,2)=0;
            d_d(i,j,2)=0;
        end
    end
    
    for i=max_x
        for j=2:1:max_y-1
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,j,1)-4*d_d(max_x-1,j,1)+d_d(max_x-2,j,1))/(2*d_x)+...
                lambda_steel*(d_d(max_x,j+1,2)-d_d(max_x,j-1,2))/(2*d_y);
            equilibrium(i,j,2)=miu_steel*(d_d(max_x,j+1,1)-d_d(max_x,j-1,1))/(2*d_y)+...
                miu_steel*(3*d_d(max_x,j,2)-4*d_d(max_x-1,j,2)+d_d(max_x-2,j,2))/(2*d_x);
        end
    end
    
    for i=2:1:max_x-1
        for j=max_y
            equilibrium(i,j,2)=lambda_steel*(d_d(i+1,max_y,1)-d_d(i-1,max_y,1))/(2*d_x)+...
                (2*miu_steel+lambda_steel)*(3*d_d(i,max_y,2)-4*d_d(i,max_y-1,2)+d_d(i,max_y-2,2))/(2*d_y);
            equilibrium(i,j,1)=miu_steel*(3*d_d(i,max_y,1)-4*d_d(i,max_y-1,1)+d_d(i,max_y-2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,max_y,2)-d_d(i-1,max_y,2))/(2*d_x);
        end
    end
    
    for i=max_x
        for j=max_y
            equilibrium(i,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
                lambda_steel*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
            equilibrium(i,j,2)=lambda_steel*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
                (2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
        end
    end
    %}
end