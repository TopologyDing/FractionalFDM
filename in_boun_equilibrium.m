function[]=in_boun_equilibrium(size_out_fic_for_in,miu_aluminum,lambda_aluminum,v_aluminum,d_x,d_y)
    global out_index_x
    global out_index_y
    global out_fic_index
    global out_dir_index
    global frac_deri_x_disp_x
    global frac_deri_x_disp_y
    global frac_deri_y_disp_x
    global frac_deri_y_disp_y
    global boun_frac_deri_x_disp_x
    global boun_frac_deri_x_disp_y
    global boun_frac_deri_y_disp_x
    global boun_frac_deri_y_disp_y
    global sec_order_frac_deri_x_disp_x
    global sec_order_frac_deri_x_disp_y
    global sec_order_frac_deri_y_disp_x
    global sec_order_frac_deri_y_disp_y
    global sym_boun_frac_deri_x_disp_x
    global sym_boun_frac_deri_x_disp_y
    global sym_boun_frac_deri_y_disp_x
    global sym_boun_frac_deri_y_disp_y
    global equilibrium
    global in_cross_boun_mark
    global d_d
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global finite_dis

    for i=2:1:size_out_fic_for_in-1
        m=out_index_x(i);
        n=out_index_y(i);
        x_local=m;
        y_local=n;
        if(in_cross_boun_mark(x_local,y_local,1)==1)
           if(in_cross_boun_mark(x_local,y_local,2)==1)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_x_disp_x(m,n)-sym_boun_frac_deri_x_disp_x(i,1))/(finite_dis(i,1))+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-frac_deri_y_disp_x(m,n-1))/(d_y)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m,n+1)-sec_order_frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
                    2*v_aluminum*(boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,1)));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_y_disp_y(m,n)-frac_deri_y_disp_y(m,n-1))/(d_y)+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-sym_boun_frac_deri_x_disp_y(i,1))/(finite_dis(i,1))+...
                    (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,1))+...
                    2*v_aluminum*(sec_order_frac_deri_x_disp_x(m,n+1)-sec_order_frac_deri_x_disp_x(m,n-1))/(2*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_x_disp_x(m,n)-frac_deri_x_disp_x(m-1,n))/(d_x)+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-sym_boun_frac_deri_y_disp_x(i,1))/(finite_dis(i,2))+...
                    (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,1)-sym_boun_frac_deri_x_disp_y(i,1))/(2*finite_dis(i,2))+...
                    2*v_aluminum*(sec_order_frac_deri_y_disp_y(m+1,n)-sec_order_frac_deri_y_disp_y(m-1,n))/(2*d_x));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_y_disp_y(m,n)-sym_boun_frac_deri_y_disp_y(i,1))/(finite_dis(i,2))+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-frac_deri_x_disp_y(m-1,n))/(d_x)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m+1,n)-sec_order_frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
                    2*v_aluminum*(boun_frac_deri_x_disp_x(i,1)-sym_boun_frac_deri_x_disp_x(i,1))/(2*finite_dis(i,2)));
           end
        elseif(in_cross_boun_mark(x_local,y_local,1)==2)
            equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_x_disp_x(m,n)-frac_deri_x_disp_x(m-1,n))/(d_x)+...
                (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-frac_deri_y_disp_x(m,n-1))/(d_y)+...
                (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m,n+1)-sec_order_frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
                2*v_aluminum*(sec_order_frac_deri_y_disp_y(m+1,n)-sec_order_frac_deri_y_disp_y(m-1,n))/(2*d_x));
            equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_y_disp_y(m,n)-frac_deri_y_disp_y(m,n-1))/(d_y)+...
                (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-frac_deri_x_disp_y(m-1,n))/(d_x)+...
                (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m+1,n)-sec_order_frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
                2*v_aluminum*(sec_order_frac_deri_x_disp_x(m,n+1)-sec_order_frac_deri_x_disp_x(m,n-1))/(2*d_y));
            
            %{
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(m-1,n,1)+d_d(m+1,n,1)-2*d_d(m,n,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n+1)-frac_deri_y_disp_x(m,n-1))/(2*d_y)+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n+1)-frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
                    2*v_aluminum*(frac_deri_y_disp_y(m+1,n)-frac_deri_y_disp_y(m-1,n))/(2*d_x));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(m,n-1,2)+d_d(m,n+1,2)-2*d_d(m,n,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m+1,n)-frac_deri_x_disp_y(m-1,n))/(2*d_x)+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m+1,n)-frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
                    2*v_aluminum*(frac_deri_x_disp_x(m,n+1)-frac_deri_x_disp_x(m,n-1))/(2*d_y));
            %}
            %{
            equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(m-1,n,1)+d_d(m+1,n,1)-2*d_d(m,n,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(m,n-1,1)+d_d(m,n+1,1)-2*d_d(m,n,1))/(d_y^2)+...
                (out_fic_for_in_uy(i,1)-d_d(m-1,n+1,2)-d_d(m+1,n-1,2)+d_d(m-1,n-1,2))/(4*d_x*d_y));
            equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(m,n-1,2)+d_d(m,n+1,2)-2*d_d(m,n,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(m-1,n,2)+d_d(m+1,n,2)-2*d_d(m,n,2))/(d_x^2)+...
                (out_fic_for_in_ux(i,1)-d_d(m-1,n+1,1)-d_d(m+1,n-1,1)+d_d(m-1,n-1,1))/(4*d_x*d_y));
            %}
            %{
            equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(m-1,n,1)+d_d(m+1,n,1)-2*d_d(m,n,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(m,n-1,1)+d_d(m,n+1,1)-2*d_d(m,n,1))/(d_y^2)+...
                (out_fic_for_in_uy(i,1)-d_d(m-1,n+1,2)-d_d(m+1,n-1,2)+d_d(m-1,n-1,2))/(4*d_x*d_y));
            equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(m,n-1,2)+d_d(m,n+1,2)-2*d_d(m,n,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(m-1,n,2)+d_d(m+1,n,2)-2*d_d(m,n,2))/(d_x^2)+...
                (out_fic_for_in_ux(i,1)-d_d(m-1,n+1,1)-d_d(m+1,n-1,1)+d_d(m-1,n-1,1))/(4*d_x*d_y));
            %}
        elseif(in_cross_boun_mark(x_local,y_local,1)==3)
            equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_x_disp_x(m,n)-sym_boun_frac_deri_x_disp_x(i,1))/(finite_dis(i,1))+...
                (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-sym_boun_frac_deri_y_disp_x(i,2))/(finite_dis(i,2))+...
                (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,2)-sym_boun_frac_deri_x_disp_y(i,2))/(2*finite_dis(i,2))+...
                2*v_aluminum*(boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,1)));
            equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_y_disp_y(m,n)-sym_boun_frac_deri_y_disp_y(i,2))/(finite_dis(i,2))+...
                (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-sym_boun_frac_deri_x_disp_y(i,1))/(finite_dis(i,1))+...
                (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,1))+...
                (2*v_aluminum)*(boun_frac_deri_x_disp_x(i,2)-sym_boun_frac_deri_x_disp_x(i,2))/(2*finite_dis(i,2)));
        elseif(in_cross_boun_mark(x_local,y_local,1)==4)
            equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_x_disp_x(m,n)-sym_boun_frac_deri_x_disp_x(i,1))/(finite_dis(i,1))+...
                (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-sym_boun_frac_deri_y_disp_x(i,2))/(finite_dis(i,2))+...
                (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,2)-sym_boun_frac_deri_x_disp_y(i,2))/(2*finite_dis(i,2))+...
                2*v_aluminum*(boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,1)));
            equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (frac_deri_y_disp_y(m,n)-sym_boun_frac_deri_y_disp_y(i,2))/(finite_dis(i,2))+...
                (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-sym_boun_frac_deri_x_disp_y(i,1))/(finite_dis(i,1))+...
                (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,1))+...
                (2*v_aluminum)*(boun_frac_deri_x_disp_x(i,2)-sym_boun_frac_deri_x_disp_x(i,2))/(2*finite_dis(i,2)));
        elseif(in_cross_boun_mark(x_local,y_local,1)==5)
            
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_x_disp_x(m,n)-frac_deri_x_disp_x(m-1,n))/(d_x)+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-sym_boun_frac_deri_y_disp_x(i,1))/(finite_dis(i,2))+...
                    (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,1)-sym_boun_frac_deri_x_disp_y(i,1))/(2*finite_dis(i,2))+...
                    2*v_aluminum*(sec_order_frac_deri_y_disp_y(m+1,n)-sec_order_frac_deri_y_disp_y(m-1,n))/(2*d_x));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_y_disp_y(m,n)-sym_boun_frac_deri_y_disp_y(i,1))/(finite_dis(i,2))+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-frac_deri_x_disp_y(m-1,n))/(d_x)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m+1,n)-sec_order_frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
                    2*v_aluminum*(boun_frac_deri_x_disp_x(i,1)-sym_boun_frac_deri_x_disp_x(i,1))/(2*finite_dis(i,2)));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_x_disp_x(m,n)-sym_boun_frac_deri_x_disp_x(i,1))/(finite_dis(i,1))+...
                    (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-frac_deri_y_disp_x(m,n-1))/(d_y)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m,n+1)-sec_order_frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
                    2*v_aluminum*(boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,1)));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (frac_deri_y_disp_y(m,n)-frac_deri_y_disp_y(m,n-1))/(d_y)+...
                    (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-sym_boun_frac_deri_x_disp_y(i,1))/(finite_dis(i,1))+...
                    (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,1))+...
                    2*v_aluminum*(sec_order_frac_deri_x_disp_x(m,n+1)-sec_order_frac_deri_x_disp_x(m,n-1))/(2*d_y));
            end
            
            %{
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (sec_order_frac_deri_x_disp_x(m+1,n)-sec_order_frac_deri_x_disp_x(m-1,n))/(2*d_x)+...
                    (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,2))+...
                    (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,1)-sym_boun_frac_deri_x_disp_y(i,1))/(2*finite_dis(i,2))+...
                    2*v_aluminum*(sec_order_frac_deri_y_disp_y(m+1,n)-sec_order_frac_deri_y_disp_y(m-1,n))/(2*d_x));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,2))+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m+1,n)-sec_order_frac_deri_x_disp_y(m-1,n))/(2*d_x)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m+1,n)-sec_order_frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
                    2*v_aluminum*(boun_frac_deri_x_disp_x(i,1)-sym_boun_frac_deri_x_disp_x(i,1))/(2*finite_dis(i,2)));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                equilibrium(m,n,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (boun_frac_deri_x_disp_x(i,1)-sym_boun_frac_deri_x_disp_x(i,1))/(2*finite_dis(i,1))+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m,n+1)-sec_order_frac_deri_y_disp_x(m,n-1))/(2*d_y)+...
                    (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m,n+1)-sec_order_frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
                    2*v_aluminum*(boun_frac_deri_y_disp_y(i,1)-sym_boun_frac_deri_y_disp_y(i,1))/(2*finite_dis(i,1)));
                equilibrium(m,n,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (sec_order_frac_deri_y_disp_y(m,n+1)-sec_order_frac_deri_y_disp_y(m,n-1))/(2*d_y)+...
                    (1-2*v_aluminum)*(boun_frac_deri_x_disp_y(i,1)-sym_boun_frac_deri_x_disp_y(i,1))/(2*finite_dis(i,1))+...
                    (1-2*v_aluminum)*(boun_frac_deri_y_disp_x(i,1)-sym_boun_frac_deri_y_disp_x(i,1))/(2*finite_dis(i,1))+...
                    2*v_aluminum*(sec_order_frac_deri_x_disp_x(m,n+1)-sec_order_frac_deri_x_disp_x(m,n-1))/(2*d_y));
            end
            %}
        end
    end
    %}
    %{
    for i=2:1:size_out_fic_for_in-1
        x_global=out_index_x(i);
        y_global=out_index_y(i);
        x_local=x_global;
        y_local=y_global;
        if(in_cross_boun_mark(x_local,y_local,1)==1)
           if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local,y_local+1);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_upper_x=out_fic_for_in_ux(j_1,1);
                    out_fic_upper_y=out_fic_for_in_uy(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==1)
                            out_fic_upper_x=out_fic_for_in_ux(j_1,k);
                            out_fic_upper_y=out_fic_for_in_uy(j_1,k);
                        end
                    end
                end
                j_2=out_fic_index(x_local,y_local-1);
                if(out_fic_for_in_ux(j_2,2)==0)
                    out_fic_bottom_x=out_fic_for_in_ux(j_2,1);
                    out_fic_bottom_y=out_fic_for_in_uy(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==1)
                            out_fic_bottom_x=out_fic_for_in_ux(j_2,k);
                            out_fic_bottom_y=out_fic_for_in_uy(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_upper_y-d_d(x_global-1,y_global+1,2)-out_fic_bottom_y+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_upper_x-d_d(x_global-1,y_global+1,1)-out_fic_bottom_x+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2) 
                j_1=out_fic_index(x_local-1,y_local);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_left_x=out_fic_for_in_ux(j_1,1);
                    out_fic_left_y=out_fic_for_in_uy(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==2)
                            out_fic_left_x=out_fic_for_in_ux(j_1,k);
                            out_fic_left_y=out_fic_for_in_uy(j_1,k);
                        end
                    end
                end
                j_2=out_fic_index(x_local+1,y_local);
                if(out_fic_for_in_ux(j_2,2)==0)
                    out_fic_right_x=out_fic_for_in_ux(j_2,1);
                    out_fic_right_y=out_fic_for_in_uy(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==2)
                            out_fic_right_x=out_fic_for_in_ux(j_2,k);
                            out_fic_right_y=out_fic_for_in_uy(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_right_y-out_fic_left_y-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_right_x-out_fic_left_x-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        elseif(in_cross_boun_mark(x_local,y_local,1)==2)
            equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (out_fic_for_in_uy(i,1)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (out_fic_for_in_ux(i,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
        elseif(in_cross_boun_mark(x_local,y_local,1)==3)
            equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global-1,y_global,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_for_in_ux(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (out_fic_for_in_uy(i,3)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global,y_global-1,2)+out_fic_for_in_uy(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (out_fic_for_in_ux(i,3)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
        elseif(in_cross_boun_mark(x_local,y_local,1)==4)
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local-1,y_local);	
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_for_in_ux(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_for_in_uy(i,3)-out_fic_for_in_uy(j_1,1)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_for_in_uy(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_for_in_ux(i,3)-out_fic_for_in_ux(j_1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                j_1=out_fic_index(x_local,y_local-1);
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_for_in_ux(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_for_in_uy(i,3)-d_d(x_global-1,y_global+1,2)-out_fic_for_in_uy(j_1,1)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_for_in_uy(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_for_in_ux(i,3)-d_d(x_global-1,y_global+1,1)-out_fic_for_in_ux(j_1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        elseif(in_cross_boun_mark(x_local,y_local,1)==5)
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local+1,y_local);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_right_x=out_fic_for_in_ux(j_1,1);
                    out_fic_right_y=out_fic_for_in_uy(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==2)
                            out_fic_right_x=out_fic_for_in_ux(j_1,k);
                            out_fic_right_y=out_fic_for_in_uy(j_1,k);
                            break;
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_right_y-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_right_x-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                j_1=out_fic_index(x_local,y_local+1);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_upper_x=out_fic_for_in_ux(j_1,1);
                    out_fic_upper_y=out_fic_for_in_uy(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==1)
                            out_fic_upper_x=out_fic_for_in_ux(j_1,k);
                            out_fic_upper_y=out_fic_for_in_uy(j_1,k);
                            break;
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_for_in_ux(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_upper_y-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_for_in_uy(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_upper_x-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        end
    end
    %}
    
    
end