function[]=in_boun_disp_grad(i,in_boun_1,in_boun_2,global_x,global_y,size_out_fic_for_in)
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global out_sur_nor_dir_x
    global out_sur_nor_dir_y
    global in_disp_x
    global in_disp_y
    global in_deri_x_disp_x
    global in_deri_x_disp_y
    global in_deri_y_disp_x
    global in_deri_y_disp_y
    global out_fic_index
    global out_dir_index
    global point_material
    global d_x
    global d_y
    global r
    global d_d
    
    center_x=0;
    center_y=0;
    if(i==1)
        dis_x=(global_x-1)*d_x;
        dis_y=sqrt(r^2-dis_x^2);
        out_sur_nor_dir_x(i,1)=dis_x;
        out_sur_nor_dir_y(i,1)=dis_y;
        %ux direction
        x_disp=[d_d(global_x,global_y-1,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,1)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
        in_deri_y_disp_x(i,1)=2*a*dis_y+b;
        %uy direction
        y_disp=[d_d(global_x,global_y-1,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,1)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
        in_deri_y_disp_y(i,1)=2*a*dis_y+b;
        
        %boundary points
        dis_x=(global_x-1)*d_x;
        dis_y=sqrt(r^2-dis_x^2); 
        %right2 situation
        j=out_fic_index(global_x+2,global_y);
        for k=1:1:2
            if(out_dir_index(j,k)==2)
                out_fic_right2_ux=out_fic_for_in_ux(j,k);
                out_fic_right2_uy=out_fic_for_in_uy(j,k);
                break;
            end
        end
        %right2 ux
        x_disp=[d_d(global_x+2,global_y-1,1);d_d(global_x+2,global_y,1);out_fic_right2_ux];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_right2_disp_x=a*dis_y^2+b*dis_y+c;
        %right2 uy
        y_disp=[d_d(global_x+2,global_y-1,2);d_d(global_x+2,global_y,2);out_fic_right2_uy];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_right2_disp_y=a*dis_y^2+b*dis_y+c;
        %right situation
        j=out_fic_index(global_x+1,global_y);
        for k=1:1:2
            if(out_dir_index(j,k)==2)
                out_fic_right_ux=out_fic_for_in_ux(j,k);
                out_fic_right_uy=out_fic_for_in_uy(j,k);
                break;
            end
        end
        %right ux
        x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);out_fic_right_ux];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_right_disp_x=a*dis_y^2+b*dis_y+c;
        %right uy
        y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);out_fic_right_uy];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_right_disp_y=a*dis_y^2+b*dis_y+c; 
        in_deri_x_disp_x(1,1)=(-3*in_disp_x(1,1)+4*in_right_disp_x-in_right2_disp_x)/(2*d_x);
        in_deri_x_disp_y(1,1)=(-3*in_disp_y(1,1)+4*in_right_disp_y-in_right2_disp_y)/(2*d_x);
    end
    if(i==length(out_fic_for_in_ux))
        dis_y=(global_y-1)*d_y;
        dis_x=sqrt(r^2-dis_y^2);
        out_sur_nor_dir_x(i,1)=dis_x;
        out_sur_nor_dir_y(i,1)=dis_y;
        %ux direction
        x_disp=[d_d(global_x-1,global_y,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,1)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
        in_deri_x_disp_x(i,1)=2*a*dis_x+b;
        %uy direction
        y_disp=[d_d(global_x-1,global_y,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,1)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
        in_deri_x_disp_y(i,1)=2*a*dis_x+b;
        
        % boundary points
        dis_y=(global_y-1)*d_y;
        dis_x=sqrt(r^2-dis_y^2);
        %upper situation
        j=out_fic_index(global_x,global_y+1);
        for k=1:1:2
            if(out_dir_index(j,k)==1)
                out_fic_upper_x=out_fic_for_in_ux(j,k);
                out_fic_upper_y=out_fic_for_in_uy(j,k);
                break;
            end
        end
        %upper ux
        x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);out_fic_upper_x];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_upper_disp_x=a*dis_x^2+b*dis_x+c;
        %upper uy            
        y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);out_fic_upper_y];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_upper_disp_y=a*dis_x^2+b*dis_x+c;
        %upper2 situation
        j=out_fic_index(global_x,global_y+2);
        for k=1:1:2
            if(out_dir_index(j,k)==1)
                out_fic_upper2_x=out_fic_for_in_ux(j,k);
                out_fic_upper2_y=out_fic_for_in_uy(j,k);
                break;
            end
        end
        %upper2 ux
        x_disp=[d_d(global_x-1,global_y+2,1);d_d(global_x,global_y+2,1);out_fic_upper2_x];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_upper2_disp_x=a*dis_x^2+b*dis_x+c;
        %upper2 uy            
        y_disp=[d_d(global_x-1,global_y+2,2);d_d(global_x,global_y+2,2);out_fic_upper2_y];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        in_upper2_disp_y=a*dis_x^2+b*dis_x+c;
        in_deri_y_disp_x(size_out_fic_for_in,1)=(-3*in_disp_x(size_out_fic_for_in,1)+4*in_upper_disp_x-in_upper2_disp_x)/(2*d_y);
        in_deri_y_disp_y(size_out_fic_for_in,1)=(-3*in_disp_y(size_out_fic_for_in,1)+4*in_upper_disp_y-in_upper2_disp_y)/(2*d_y);
    end
    if((in_boun_1==1)||(in_boun_1==5))
        if(((in_boun_1==1)&&(in_boun_2==1))||((in_boun_1==5)&&((in_boun_2==2))))
            dis_y=(global_y-1)*d_y-center_y;
            dis_x=sqrt(r^2-dis_y^2)+center_x;
            out_sur_nor_dir_x(i,1)=dis_x-center_x;
            out_sur_nor_dir_y(i,1)=dis_y;
            %ux direction
            x_disp=[d_d(global_x-1,global_y,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
            in_deri_x_disp_x(i,1)=2*a*dis_x+b;
            %uy direction
            y_disp=[d_d(global_x-1,global_y,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
            in_deri_x_disp_y(i,1)=2*a*dis_x+b;
        elseif(((in_boun_1==1)&&(in_boun_2==2))||((in_boun_1==5)&&((in_boun_2==1))))
            dis_x=(global_x-1)*d_x-center_x;
            dis_y=sqrt(r^2-dis_x^2)+center_y;
            out_sur_nor_dir_x(i,1)=dis_x;
            out_sur_nor_dir_y(i,1)=dis_y-center_y;
            %ux direction
            x_disp=[d_d(global_x,global_y-1,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
            in_deri_y_disp_x(i,1)=2*a*dis_y+b;
            %uy direction
            y_disp=[d_d(global_x,global_y-1,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
            in_deri_y_disp_y(i,1)=2*a*dis_y+b;
        end
    elseif((in_boun_1==3)||(in_boun_1==4))
        if(((in_boun_1==3)&&(in_boun_2==1))||((in_boun_1==4)&&((in_boun_2==1)||(in_boun_2==2))))
            dis_y=(global_y-1)*d_y-center_y;
            dis_x=sqrt(r^2-dis_y^2)+center_x;
            out_sur_nor_dir_x(i,1)=dis_x-center_x;
            out_sur_nor_dir_y(i,1)=dis_y;
            %ux direction
            x_disp=[d_d(global_x-1,global_y,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
            in_deri_x_disp_x(i,1)=2*a*dis_x+b;
            %uy direction
            y_disp=[d_d(global_x-1,global_y,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
            in_deri_x_disp_y(i,1)=2*a*dis_x+b; 
            
            dis_x=(global_x-1)*d_x-center_x;
            dis_y=sqrt(r^2-dis_x^2)+center_y;
            out_sur_nor_dir_x(i,2)=dis_x;
            out_sur_nor_dir_y(i,2)=dis_y-center_y;
            %ux direction
            x_disp=[d_d(global_x,global_y-1,1);d_d(global_x,global_y,1);out_fic_for_in_ux(i,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_x(i,2)=a*dis_y^2+b*dis_y+c;
            in_deri_y_disp_x(i,2)=2*a*dis_y+b;
            %uy direction
            y_disp=[d_d(global_x,global_y-1,2);d_d(global_x,global_y,2);out_fic_for_in_uy(i,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_disp_y(i,2)=a*dis_y^2+b*dis_y+c;
            in_deri_y_disp_y(i,2)=2*a*dis_y+b;
        end
    elseif(in_boun_1==2)
    end 

    if(in_boun_1==1)
        if(in_boun_2==1)
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper situation
            j=out_fic_index(global_x,global_y+1);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_upper_x=out_fic_for_in_ux(j,1);
                out_fic_upper_y=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==1)
                        out_fic_upper_x=out_fic_for_in_ux(j,k);
                        out_fic_upper_y=out_fic_for_in_uy(j,k);
                    end
                end
            end
            %upper ux
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);out_fic_upper_x];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %upper uy            
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);out_fic_upper_y];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_y=a*dis_x^2+b*dis_x+c;
            %bottom situation
            j=out_fic_index(global_x,global_y-1);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_bottom_x=out_fic_for_in_ux(j,1);
                out_fic_bottom_y=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==1)
                        out_fic_bottom_x=out_fic_for_in_ux(j,k);
                        out_fic_bottom_y=out_fic_for_in_uy(j,k);
                    end
                end
            end
            %bottom ux
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);out_fic_bottom_x];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;   
            %bottom uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);out_fic_bottom_y];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;  
            in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
        elseif(in_boun_2==2)
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2); 
            %left situation
            j=out_fic_index(global_x-1,global_y);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_left_ux=out_fic_for_in_ux(j,1);
                out_fic_left_uy=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==2)
                        out_fic_left_ux=out_fic_for_in_ux(j,k);
                        out_fic_left_uy=out_fic_for_in_uy(j,k);
                    end
                end
            end
            %left ux
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);out_fic_left_ux];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            %left uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);out_fic_left_uy];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            %right situation
            j=out_fic_index(global_x+1,global_y);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_right_ux=out_fic_for_in_ux(j,1);
                out_fic_right_uy=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==2)
                        out_fic_right_ux=out_fic_for_in_ux(j,k);
                        out_fic_right_uy=out_fic_for_in_uy(j,k);
                    end
                end
            end
            %right ux
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);out_fic_right_ux];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_x=a*dis_y^2+b*dis_y+c;
            %right uy
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);out_fic_right_uy];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_y=a*dis_y^2+b*dis_y+c; 
            in_deri_x_disp_x(i,1)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,1)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
        end
    elseif(in_boun_1==5)
        if(in_boun_2==1)
            %left2 situation
            %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            x_disp=[d_d(global_x-2,global_y-1,1);d_d(global_x-2,global_y,1);d_d(global_x-2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x-2,global_y-1,2);d_d(global_x-2,global_y,2);d_d(global_x-2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_y=a*dis_y^2+b*dis_y+c;
            %left situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
            %}
            
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            j=out_fic_index(global_x+1,global_y);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_right_x=out_fic_for_in_ux(j,1);
                out_fic_right_y=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==2)
                        out_fic_right_x=out_fic_for_in_ux(j,k);
                        out_fic_right_y=out_fic_for_in_uy(j,k);
                    end
                end
            end
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);out_fic_right_x];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);out_fic_right_y];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_y=a*dis_y^2+b*dis_y+c;
            %left situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,1)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,1)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
            %}
        elseif(in_boun_2==2)
            %bottom2 situation
            
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            j=out_fic_index(global_x,global_y+1);
            if(out_fic_for_in_ux(j,2)==0)
                out_fic_upper_x=out_fic_for_in_ux(j,1);
                out_fic_upper_y=out_fic_for_in_uy(j,1);
            elseif(out_fic_for_in_ux(j,2)~=0)
                for k=1:1:2
                    if(out_dir_index(j,k)==1)
                        out_fic_upper_x=out_fic_for_in_ux(j,k);
                        out_fic_upper_y=out_fic_for_in_uy(j,k);
                    end
                end
            end
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);out_fic_upper_x];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);out_fic_upper_y];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_y=a*dis_x^2+b*dis_x+c;  
            %bottom situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
            %}
            %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            x_disp=[d_d(global_x-1,global_y-2,1);d_d(global_x,global_y-2,1);d_d(global_x+1,global_y-2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y-2,2);d_d(global_x,global_y-2,2);d_d(global_x+1,global_y-2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_y=a*dis_x^2+b*dis_x+c; 
            %bottom situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
            %}
        end
    elseif(in_boun_1==3)
        %right situation:('R')
        %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %bottom situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            %bottom2 situation
            x_disp=[d_d(global_x-1,global_y-2,1);d_d(global_x,global_y-2,1);d_d(global_x+1,global_y-2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-2,2);d_d(global_x,global_y-2,2);d_d(global_x+1,global_y-2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
        %upper situation:('T')
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %left situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;  
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;  
            %left2 situation
            x_disp=[d_d(global_x-2,global_y-1,1);d_d(global_x-2,global_y,1);d_d(global_x-2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_x=a*dis_y^2+b*dis_y+c;  
            y_disp=[d_d(global_x-2,global_y-1,2);d_d(global_x-2,global_y,2);d_d(global_x-2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
            %}
            
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %bottom situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            %top situation
            x_disp=[d_d(global_x-1,global_y+1,1);out_fic_for_in_ux(i,2);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y+1,2);out_fic_for_in_uy(i,2);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
        %upper situation:('T')
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %left situation
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;  
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;  
            %right situation
            x_disp=[d_d(global_x+1,global_y-1,1);out_fic_for_in_ux(i,1);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_x=a*dis_y^2+b*dis_y+c;  
            y_disp=[d_d(global_x+1,global_y-1,2);out_fic_for_in_uy(i,1);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,2)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,2)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
            %}
    elseif(in_boun_1==4)
        if(in_boun_2==1)
        %right situation:('R')
        %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %bottom
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;   
            %bottom2
            x_disp=[d_d(global_x-1,global_y-2,1);d_d(global_x,global_y-2,1);d_d(global_x+1,global_y-2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-2,2);d_d(global_x,global_y-2,2);d_d(global_x+1,global_y-2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
            %}
            
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %bottom
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;   
            %upper
            j=out_fic_index(global_x-1,global_y);
            x_disp=[out_fic_for_in_ux(j,1);out_fic_for_in_ux(i,2);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[out_fic_for_in_uy(j,1);out_fic_for_in_uy(i,2);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
            %}
        %upper situation:('T')
        %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2); 
            %left
            j_1=out_fic_index(global_x-1,global_y);
            j_2=out_fic_index(global_x-2,global_y);
            if(point_material(global_x-2,global_y+1)==1)
                x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);out_fic_for_in_ux(j_1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);out_fic_for_in_uy(j_1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_y=a*dis_y^2+b*dis_y+c;
                if(out_fic_for_in_ux(j_2,2)==0)
                    out_left2_x=out_fic_for_in_ux(j_2,1);
                    out_left2_y=out_fic_for_in_uy(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==2)
                            out_left2_x=out_fic_for_in_ux(j_2,k);
                            out_left2_y=out_fic_for_in_uy(j_2,k);
                        end
                    end
                end
                x_disp=[d_d(global_x-2,global_y-1,1);d_d(global_x-2,global_y,1);out_left2_x];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left2_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[d_d(global_x-2,global_y-1,2);d_d(global_x-2,global_y,2);out_left2_y];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left2_disp_y=a*dis_y^2+b*dis_y+c;
                in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
                in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
            elseif(point_material(global_x-2,global_y+1)==-1)
				for k=1:1:2
					if(out_dir_index(j_1,k)==2)
						out_left_x=out_fic_for_in_ux(j_1,k);
						out_left_y=out_fic_for_in_uy(j_1,k);
					end
				end
                x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);out_left_x];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);out_left_y];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_y=a*dis_y^2+b*dis_y+c;
                x_disp=[d_d(global_x-2,global_y-1,1);d_d(global_x-2,global_y,1);d_d(global_x-2,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left2_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[d_d(global_x-2,global_y-1,2);d_d(global_x-2,global_y,2);d_d(global_x-2,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left2_disp_y=a*dis_y^2+b*dis_y+c;
                in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
                in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
            end
            %}
        
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2); 
            x_disp=[d_d(global_x+1,global_y-1,1);out_fic_for_in_ux(i,1);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x+1,global_y-1,2);out_fic_for_in_uy(i,1);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_y=a*dis_y^2+b*dis_y+c;
            j=out_fic_index(global_x-1,global_y);
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);out_fic_for_in_ux(j,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);out_fic_for_in_uy(j,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,2)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,2)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
            %}
        elseif(in_boun_2==2)
        %right situation:('R')
        %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %bottom
            j_1=out_fic_index(global_x,global_y-1);
            j_2=out_fic_index(global_x,global_y-2);
            if(point_material(global_x+1,global_y-2)==1)
                x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);out_fic_for_in_ux(j_1,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);out_fic_for_in_uy(j_1,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                if(out_fic_for_in_ux(j_2,2)==0)
					out_bottom2_x=out_fic_for_in_ux(j_2,1);
					out_bottom2_y=out_fic_for_in_uy(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==1)
							out_bottom2_x=out_fic_for_in_ux(j_2,k);
							out_bottom2_y=out_fic_for_in_uy(j_2,k);
                        end
                    end
                end
                x_disp=[d_d(global_x-1,global_y-2,1);d_d(global_x,global_y-2,1);out_bottom2_x];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(global_x-1,global_y-2,2);d_d(global_x,global_y-2,2);out_bottom2_y];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
                in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
                in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
            elseif(point_material(global_x+1,global_y-2)==-1)
				for k=1:1:2
					if(out_dir_index(j_1,k)==1)
						out_bottom_x=out_fic_for_in_ux(j_1,k);
						out_bottom_y=out_fic_for_in_uy(j_1,k);
					end
				end
                x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);out_bottom_x];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);out_bottom_y];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                x_disp=[d_d(global_x-1,global_y-2,1);d_d(global_x,global_y-2,1);d_d(global_x+1,global_y-2,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(global_x-1,global_y-2,2);d_d(global_x,global_y-2,2);d_d(global_x+1,global_y-2,2)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
                in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
                in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
            end
        %}
        
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            j=out_fic_index(global_x,global_y-1);
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x,global_y-1,1);out_fic_for_in_ux(j,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x,global_y-1,2);out_fic_for_in_uy(j,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            x_disp=[d_d(global_x-1,global_y+1,1);out_fic_for_in_ux(i,2);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[d_d(global_x-1,global_y+1,2);out_fic_for_in_uy(i,2);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_upper_disp_y=a*dis_x^2+b*dis_x+c;
            in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
            in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
        %}
        %upper situation:('T') 
        %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2); 
            %left
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            %left2
            x_disp=[d_d(global_x-2,global_y-1,1);d_d(global_x-2,global_y,1);d_d(global_x-2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x-2,global_y-1,2);d_d(global_x-2,global_y,2);d_d(global_x-2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left2_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
            %}
            
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2); 
            %left
            x_disp=[d_d(global_x-1,global_y-1,1);d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[d_d(global_x-1,global_y-1,2);d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_left_disp_y=a*dis_y^2+b*dis_y+c;
            %right
            j=out_fic_index(global_x,global_y-1);
            x_disp=[out_fic_for_in_ux(j,1);out_fic_for_in_ux(i,1);out_fic_for_in_ux(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[out_fic_for_in_uy(j,1);out_fic_for_in_uy(i,1);out_fic_for_in_uy(i,3)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_right_disp_y=a*dis_y^2+b*dis_y+c;
            in_deri_x_disp_x(i,2)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
            in_deri_x_disp_y(i,2)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
            %}
        end
    end
end