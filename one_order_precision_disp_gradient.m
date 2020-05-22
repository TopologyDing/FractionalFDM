function[]=one_order_precision_disp_gradient(high_x,high_y,d_x,d_y)
    global one_order_disp_gradient_x
    global one_order_disp_gradient_y
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global out_fic_index
    global in_cross_boun_mark
    global d_d
    global point_material
    
    for i=1:1:high_x
        for j=1:1:high_y
            if(point_material(i,j)==-2)
                one_order_disp_gradient_x(i,j,1)=(d_d(i+1,j,1)-d_d(i,j,1))/d_x;
                one_order_disp_gradient_x(i,j,2)=(d_d(i+1,j,2)-d_d(i,j,2))/d_x;
                one_order_disp_gradient_y(i,j,1)=(d_d(i,j+1,1)-d_d(i,j,1))/d_y;
                one_order_disp_gradient_y(i,j,2)=(d_d(i,j+1,2)-d_d(i,j,2))/d_y;
            elseif(point_material(i,j)==-1)
                t=out_fic_index(i,j);
                boun_1=in_cross_boun_mark(i,j,1);
                boun_2=in_cross_boun_mark(i,j,2);
                if(boun_1==-1)
                    if(boun_2==1)
                        one_order_disp_gradient_x(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(d_d(i,j+1,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(d_d(i,j+1,2)-d_d(i,j,2))/d_y;
                    else
                        one_order_disp_gradient_x(i,j,1)=(d_d(i+1,j,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(d_d(i+1,j,2)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_y;
                    end
                elseif(boun_1==1)
                    if(boun_2==1)
                        one_order_disp_gradient_x(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(d_d(i,j+1,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(d_d(i,j+1,2)-d_d(i,j,2))/d_y;
                    else
                        one_order_disp_gradient_x(i,j,1)=(d_d(i+1,j,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(d_d(i+1,j,2)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_y;
                    end
                elseif(boun_1==2)
                    one_order_disp_gradient_x(i,j,1)=(d_d(i+1,j,1)-d_d(i,j,1))/d_x;
                    one_order_disp_gradient_x(i,j,2)=(d_d(i+1,j,2)-d_d(i,j,2))/d_x;
                    one_order_disp_gradient_y(i,j,1)=(d_d(i,j+1,1)-d_d(i,j,1))/d_y;
                    one_order_disp_gradient_y(i,j,2)=(d_d(i,j+1,2)-d_d(i,j,2))/d_y;
                elseif(boun_1==3)
                    one_order_disp_gradient_x(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_x;
                    one_order_disp_gradient_x(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_x;
                    one_order_disp_gradient_y(i,j,1)=(out_fic_for_in_ux(t,2)-d_d(i,j,1))/d_y;
                    one_order_disp_gradient_y(i,j,2)=(out_fic_for_in_uy(t,2)-d_d(i,j,2))/d_y;
                elseif(boun_1==4)
                    one_order_disp_gradient_x(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_x;
                    one_order_disp_gradient_x(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_x;
                    one_order_disp_gradient_y(i,j,1)=(out_fic_for_in_ux(t,2)-d_d(i,j,1))/d_y;
                    one_order_disp_gradient_y(i,j,2)=(out_fic_for_in_uy(t,2)-d_d(i,j,2))/d_y;
                elseif(boun_1==5)
                    if(boun_2==2)
                        one_order_disp_gradient_x(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(d_d(i,j+1,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(d_d(i,j+1,2)-d_d(i,j,2))/d_y;
                    else
                        one_order_disp_gradient_x(i,j,1)=(d_d(i+1,j,1)-d_d(i,j,1))/d_x;
                        one_order_disp_gradient_x(i,j,2)=(d_d(i+1,j,2)-d_d(i,j,2))/d_x;
                        one_order_disp_gradient_y(i,j,1)=(out_fic_for_in_ux(t,1)-d_d(i,j,1))/d_y;
                        one_order_disp_gradient_y(i,j,2)=(out_fic_for_in_uy(t,1)-d_d(i,j,2))/d_y;
                    end
                end
            end
        end
    end
end