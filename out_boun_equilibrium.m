function[]=out_boun_equilibrium(size_in_fic_for_out,miu_steel,lambda_steel,v_steel,d_x,d_y)
    global out_cross_boun_mark
    global in_fic_for_out_ux
    global in_fic_for_out_uy
    global in_index_x
    global in_index_y
    global in_dir_index
    global in_fic_index
    global equilibrium
    global d_d
    for i=2:1:size_in_fic_for_out-1
        x_global=in_index_x(i);
        y_global=in_index_y(i);
        x_local=x_global;
        y_local=y_global;
        if(out_cross_boun_mark(x_local,y_local,1)==1)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local,y_local+1);
                if(in_fic_for_out_ux(j_1,2)==0)
                    in_fic_upper_x=in_fic_for_out_ux(j_1,1);
                    in_fic_upper_y=in_fic_for_out_uy(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==1)
                            in_fic_upper_x=in_fic_for_out_ux(j_1,k);
                            in_fic_upper_y=in_fic_for_out_uy(j_1,k);
                        end
                    end
                end
                j_2=in_fic_index(x_local,y_local-1);
                if(in_fic_for_out_ux(j_2,2)==0)
                    in_fic_bottom_x=in_fic_for_out_ux(j_2,1);
                    in_fic_bottom_y=in_fic_for_out_uy(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==1)
                            in_fic_bottom_x=in_fic_for_out_ux(j_2,k);
                            in_fic_bottom_y=in_fic_for_out_uy(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_ux(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-in_fic_upper_y-d_d(x_global+1,y_global-1,2)+in_fic_bottom_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_for_out_uy(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-in_fic_upper_x-d_d(x_global+1,y_global-1,1)+in_fic_bottom_x)/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local-1,y_local);
                if(in_fic_for_out_ux(j_1,2)==0)
                    in_fic_left_x=in_fic_for_out_ux(j_1,1);
                    in_fic_left_y=in_fic_for_out_uy(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==2)
                            in_fic_left_x=in_fic_for_out_ux(j_1,k);
                            in_fic_left_y=in_fic_for_out_uy(j_1,k);
                        end
                    end
                end
                j_2=in_fic_index(x_local+1,y_local);
                if(in_fic_for_out_ux(j_2,2)==0)
                    in_fic_right_x=in_fic_for_out_ux(j_2,1);
                    in_fic_right_y=in_fic_for_out_uy(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==2)
                            in_fic_right_x=in_fic_for_out_ux(j_2,k);
                            in_fic_right_y=in_fic_for_out_uy(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_for_out_ux(i,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-in_fic_right_y+in_fic_left_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_uy(i,1)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-in_fic_right_x+in_fic_left_x)/(4*d_x*d_y));
            end
        elseif(out_cross_boun_mark(x_local,y_local,1)==2)
            equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (in_fic_for_out_ux(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_steel)*(in_fic_for_out_ux(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_for_out_uy(i,3))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (in_fic_for_out_uy(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_steel)*(in_fic_for_out_uy(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_for_out_ux(i,3))/(4*d_x*d_y));
        elseif(out_cross_boun_mark(x_local,y_local,1)==3)
            equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_for_out_uy(i,1))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_for_out_ux(i,1))/(4*d_x*d_y));
        elseif(out_cross_boun_mark(x_local,y_local,1)==4)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local-1,y_local);
                if(in_fic_for_out_ux(j_1,2)==0)
                    in_fic_left_x=in_fic_for_out_ux(j_1,1);
                    in_fic_left_y=in_fic_for_out_uy(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==2)
                            in_fic_left_x=in_fic_for_out_ux(j_1,k);
                            in_fic_left_y=in_fic_for_out_uy(j_1,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_for_out_ux(i,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_left_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_uy(i,1)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_left_x)/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local,y_local-1);
                if(in_fic_for_out_ux(j_1,2)==0)
                    in_fic_bottom_x=in_fic_for_out_ux(j_1,1);
                    in_fic_bottom_y=in_fic_for_out_uy(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==1)
                            in_fic_bottom_x=in_fic_for_out_ux(j_1,k);
                            in_fic_bottom_y=in_fic_for_out_uy(j_1,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_ux(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_bottom_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_for_out_uy(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_bottom_x)/(4*d_x*d_y));
            end
        elseif(out_cross_boun_mark(x_local,y_local,1)==5)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local+1,y_local);
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_ux(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_for_out_ux(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-in_fic_for_out_uy(j_1,1)+in_fic_for_out_uy(i,3))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_uy(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_for_out_uy(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-in_fic_for_out_ux(j_1,1)+in_fic_for_out_ux(i,3))/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local,y_local+1);
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_ux(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_for_out_ux(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-in_fic_for_out_uy(j_1,1)-d_d(x_global+1,y_global-1,2)+in_fic_for_out_uy(i,3))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_for_out_uy(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_for_out_uy(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-in_fic_for_out_ux(j_1,1)-d_d(x_global+1,y_global-1,1)+in_fic_for_out_ux(i,3))/(4*d_x*d_y));
            end

        end
    end
end