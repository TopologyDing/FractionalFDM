function[]=normal_point_equilibrium(max_x,max_y,miu_aluminum,lambda_aluminum,v_aluminum,miu_steel,lambda_steel,v_steel)
    global equilibrium
    global point_material
    global frac_deri_x_disp_x
    global frac_deri_x_disp_y
    global frac_deri_y_disp_x
    global frac_deri_y_disp_y
    global sec_order_frac_deri_x_disp_x
    global sec_order_frac_deri_x_disp_y
    global sec_order_frac_deri_y_disp_x
    global sec_order_frac_deri_y_disp_y
    global d_d
    global d_x
    global d_y
    
    for i=2:1:max_x-1
        for j=2:1:max_y-1
            if(point_material(i,j)==-2)
                equilibrium(i,j,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_x_disp_x(i,j)-frac_deri_x_disp_x(i-1,j))/(d_x)+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(i,j)-frac_deri_y_disp_x(i,j-1))/(d_y)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(i,j+1)-sec_order_frac_deri_x_disp_y(i,j-1))/(2*d_y)+...
                    2*v_aluminum*(sec_order_frac_deri_y_disp_y(i+1,j)-sec_order_frac_deri_y_disp_y(i-1,j))/(2*d_x));
                equilibrium(i,j,2)=(miu_aluminum+lambda_aluminum)*((1-2*v_aluminum)*...
                    (frac_deri_x_disp_y(i,j)-frac_deri_x_disp_y(i-1,j))/(d_x)+...
                    2*(1-v_aluminum)*(frac_deri_y_disp_y(i,j)-frac_deri_y_disp_y(i,j-1))/(d_y)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(i+1,j)-sec_order_frac_deri_y_disp_x(i-1,j))/(2*d_x)+...
                    2*v_aluminum*(sec_order_frac_deri_x_disp_x(i,j+1)-sec_order_frac_deri_x_disp_x(i,j-1))/(2*d_y));
            elseif(point_material(i,j)==2)
                equilibrium(i,j,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*(d_d(i-1,j,1)+d_d(i+1,j,1)-2*d_d(i,j,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(i,j-1,1)+d_d(i,j+1,1)-2*d_d(i,j,1))/(d_y^2)+...
                    (d_d(i+1,j+1,2)-d_d(i-1,j+1,2)-d_d(i+1,j-1,2)+d_d(i-1,j-1,2))/(4*d_x*d_y));
                equilibrium(i,j,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*(d_d(i,j-1,2)+d_d(i,j+1,2)-2*d_d(i,j,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(i-1,j,2)+d_d(i+1,j,2)-2*d_d(i,j,2))/(d_x^2)+...
                    (d_d(i+1,j+1,1)-d_d(i-1,j+1,1)-d_d(i+1,j-1,1)+d_d(i-1,j-1,1))/(4*d_x*d_y));
            end
        end
    end
    %}
    %{
    for i=2:1:max_x-1
        for j=2:1:max_y-1
            if(point_material(i,j)==-2)
                equilibrium(i,j,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*(d_d(i-1,j,1)+d_d(i+1,j,1)-2*d_d(i,j,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(i,j-1,1)+d_d(i,j+1,1)-2*d_d(i,j,1))/(d_y^2)+...
                    (d_d(i+1,j+1,2)-d_d(i-1,j+1,2)-d_d(i+1,j-1,2)+d_d(i-1,j-1,2))/(4*d_x*d_y));
                equilibrium(i,j,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*(d_d(i,j-1,2)+d_d(i,j+1,2)-2*d_d(i,j,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(i-1,j,2)+d_d(i+1,j,2)-2*d_d(i,j,2))/(d_x^2)+...
                    (d_d(i+1,j+1,1)-d_d(i-1,j+1,1)-d_d(i+1,j-1,1)+d_d(i-1,j-1,1))/(4*d_x*d_y));
            elseif(point_material(i,j)==2)
                equilibrium(i,j,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*(d_d(i-1,j,1)+d_d(i+1,j,1)-2*d_d(i,j,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(i,j-1,1)+d_d(i,j+1,1)-2*d_d(i,j,1))/(d_y^2)+...
                    (d_d(i+1,j+1,2)-d_d(i-1,j+1,2)-d_d(i+1,j-1,2)+d_d(i-1,j-1,2))/(4*d_x*d_y));
                equilibrium(i,j,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*(d_d(i,j-1,2)+d_d(i,j+1,2)-2*d_d(i,j,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(i-1,j,2)+d_d(i+1,j,2)-2*d_d(i,j,2))/(d_x^2)+...
                    (d_d(i+1,j+1,1)-d_d(i-1,j+1,1)-d_d(i+1,j-1,1)+d_d(i-1,j-1,1))/(4*d_x*d_y));
            end
        end
    end
    %}
end