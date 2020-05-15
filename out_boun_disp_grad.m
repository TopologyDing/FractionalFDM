function[]=out_boun_disp_grad(i,out_boun_1,out_boun_2,global_x,global_y,size_in_fic_for_out)
    global in_fic_for_out_ux
    global in_fic_for_out_uy
    global out_disp_x
    global out_disp_y
    global out_deri_x_disp_x
    global out_deri_x_disp_y
    global out_deri_y_disp_x
    global out_deri_y_disp_y
    global in_fic_index
    global in_dir_index
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
        x_disp=[in_fic_for_out_ux(i,1);d_d(global_x,global_y,1);d_d(global_x,global_y+1,1)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
        out_deri_y_disp_x(i,1)=2*a*dis_y+b;
        %uy situation
        y_disp=[in_fic_for_out_uy(i,1);d_d(global_x,global_y,2);d_d(global_x,global_y+1,2)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
        out_deri_y_disp_y(i,1)=2*a*dis_y+b;
    end
    if(i==length(in_fic_for_out_ux))
        dis_y=(global_y-1)*d_y;
        dis_x=sqrt(r^2-dis_y^2);
        %ux situation
        x_disp=[in_fic_for_out_ux(i,1);d_d(global_x,global_y,1);d_d(global_x+1,global_y,1)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
        out_deri_x_disp_x(i,1)=2*a*dis_x+b;
        %uy situation
        y_disp=[in_fic_for_out_uy(i,1);d_d(global_x,global_y,2);d_d(global_x+1,global_y,2)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
        out_deri_x_disp_y(i,1)=2*a*dis_x+b; 
    end
    if((out_boun_1==1)||(out_boun_1==4))
        if(((out_boun_1==1)&&(out_boun_2==1))||((out_boun_1==4)&&((out_boun_2==2))))
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %ux situation
            x_disp=[in_fic_for_out_ux(i,1);d_d(global_x,global_y,1);d_d(global_x+1,global_y,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
            out_deri_x_disp_x(i,1)=2*a*dis_x+b;
            %uy situation
            y_disp=[in_fic_for_out_uy(i,1);d_d(global_x,global_y,2);d_d(global_x+1,global_y,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
            out_deri_x_disp_y(i,1)=2*a*dis_x+b; 
        elseif(((out_boun_1==1)&&(out_boun_2==2))||((out_boun_1==4)&&((out_boun_2==1))))
            dis_x=(global_x-1)*d_x-center_x;
            dis_y=sqrt(r^2-dis_x^2)+center_y;
            x_disp=[in_fic_for_out_ux(i,1);d_d(global_x,global_y,1);d_d(global_x,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
            out_deri_y_disp_x(i,1)=2*a*dis_y+b;
            %uy situation
            y_disp=[in_fic_for_out_uy(i,1);d_d(global_x,global_y,2);d_d(global_x,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
            out_deri_y_disp_y(i,1)=2*a*dis_y+b;
        end
    elseif((out_boun_1==2)||(out_boun_1==5))
        if(((out_boun_1==2)&&(out_boun_2==1))||((out_boun_1==5)&&((out_boun_2==1)||(out_boun_2==2))))
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %ux situation
            x_disp=[in_fic_for_out_ux(i,1);d_d(global_x,global_y,1);d_d(global_x+1,global_y,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
            out_deri_x_disp_x(i,1)=2*a*dis_x+b;
            %uy situation
            y_disp=[in_fic_for_out_uy(i,1);d_d(global_x,global_y,2);d_d(global_x+1,global_y,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
            out_deri_x_disp_y(i,1)=2*a*dis_x+b; 

            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            x_disp=[in_fic_for_out_ux(i,2);d_d(global_x,global_y,1);d_d(global_x,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_x(i,2)=a*dis_y^2+b*dis_y+c;
            out_deri_y_disp_x(i,2)=2*a*dis_y+b;
            %uy situation
            y_disp=[in_fic_for_out_uy(i,2);d_d(global_x,global_y,2);d_d(global_x,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_disp_y(i,2)=a*dis_y^2+b*dis_y+c;
            out_deri_y_disp_y(i,2)=2*a*dis_y+b;
        end
    elseif(out_boun_1==3)
    end
    if(i==1)
        dis_x=(global_x-1)*d_x;
        dis_y=sqrt(r^2-dis_x^2);
        %right2 situation
        j=in_fic_index(global_x+2,global_y);
        if(in_fic_for_out_ux(j,2)==0)
            in_fic_right2_x=in_fic_for_out_ux(j,1);
            in_fic_right2_y=in_fic_for_out_uy(j,1);
        elseif(in_fic_for_out_ux(j,2)~=0)
            for k=1:1:2
                if(in_dir_index(j,k)==2)
                    in_fic_right2_x=in_fic_for_out_ux(j,k);
                    in_fic_right2_y=in_fic_for_out_uy(j,k);
                end
            end
        end
        %compute right2 out_disp_x
        x_disp=[in_fic_right2_x;d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_right2_disp_x=a*dis_y^2+b*dis_y+c;
        %compute right2 out_disp_y
        y_disp=[in_fic_right2_y;d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_right2_disp_y=a*dis_y^2+b*dis_y+c;
        %right situation
        j=in_fic_index(global_x+1,global_y);
        if(in_fic_for_out_ux(j,2)==0)
            in_fic_right_x=in_fic_for_out_ux(j,1);
            in_fic_right_y=in_fic_for_out_uy(j,1);
        elseif(in_fic_for_out_ux(j,2)~=0)
            for k=1:1:2
                if(in_dir_index(j,k)==2)
                    in_fic_right_x=in_fic_for_out_ux(j,k);
                    in_fic_right_y=in_fic_for_out_uy(j,k);
                end
            end
        end
        %compute right out_disp_x
        x_disp=[in_fic_right_x;d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_right_disp_x=a*dis_y^2+b*dis_y+c;
        %compute right out_disp_y
        y_disp=[in_fic_right_y;d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
        syms a b c;
        eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
        eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
        eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_right_disp_y=a*dis_y^2+b*dis_y+c;
        out_deri_x_disp_x(1,1)=(-3*out_disp_x(1,1)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
        out_deri_x_disp_y(1,1)=(-3*out_disp_y(1,1)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
    elseif(i==size_in_fic_for_out)
        dis_y=(global_y-1)*d_y;
        dis_x=sqrt(r^2-dis_y^2);
        %the index of in_fic on outer point (x,y+1)
        %the corresponding in_fic_point is on the left of (x,y+1)
        %upper situation
        j=in_fic_index(global_x,global_y+1);
        if(in_fic_for_out_ux(j,2)==0)
            in_fic_upper_x=in_fic_for_out_ux(j,1);
            in_fic_upper_y=in_fic_for_out_uy(j,1);
        elseif(in_fic_for_out_ux(j,2)~=0)
            for k=1:1:2
                if(in_dir_index(j,k)==1)
                    in_fic_upper_x=in_fic_for_out_ux(j,k);
                    in_fic_upper_y=in_fic_for_out_uy(j,k);
                end
            end
        end
        %compute upper out_disp_x
        x_disp=[in_fic_upper_x;d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_upper_disp_x=a*dis_x^2+b*dis_x+c;
        %compute upper out_disp_y
        y_disp=[in_fic_upper_y;d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_upper_disp_y=a*dis_x^2+b*dis_x+c;
        % upper2 situation
        j=in_fic_index(global_x,global_y+2);
        if(in_fic_for_out_ux(j,2)==0)
            in_fic_upper2_x=in_fic_for_out_ux(j,1);
            in_fic_upper2_y=in_fic_for_out_uy(j,1);
        elseif(in_fic_for_out_ux(j,2)~=0)
            for k=1:1:2
                if(in_dir_index(j,k)==1)
                    in_fic_upper2_x=in_fic_for_out_ux(j,k);
                    in_fic_upper2_y=in_fic_for_out_uy(j,k);
                end
            end
        end
        %compute upper2 out_disp_x
        x_disp=[in_fic_upper2_x;d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
        %compute upper2 out_disp_y            
        y_disp=[in_fic_upper2_y;d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
        syms a b c;
        eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
        eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
        eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        out_upper2_disp_y=a*dis_x^2+b*dis_x+c;  
        %compute y direction derivative
        out_deri_y_disp_x(size_in_fic_for_out,1)=(-3*out_disp_x(size_in_fic_for_out,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
        out_deri_y_disp_y(size_in_fic_for_out,1)=(-3*out_disp_y(size_in_fic_for_out,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
    end
    if(out_boun_1==1)
        if(out_boun_2==1)
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %the index of in_fic on outer point (x,y+1)
            %the corresponding in_fic_point is on the left of (x,y+1)
            %upper situation
            j=in_fic_index(global_x,global_y+1);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_upper_x=in_fic_for_out_ux(j,1);
                in_fic_upper_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==1)
                        in_fic_upper_x=in_fic_for_out_ux(j,k);
                        in_fic_upper_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute upper out_disp_x
            x_disp=[in_fic_upper_x;d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %compute upper out_disp_y
            y_disp=[in_fic_upper_y;d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            % bottom situation
            j=in_fic_index(global_x,global_y-1);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_bottom_x=in_fic_for_out_ux(j,1);
                in_fic_bottom_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==1)
                        in_fic_bottom_x=in_fic_for_out_ux(j,k);
                        in_fic_bottom_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute bottom out_disp_x
            x_disp=[in_fic_bottom_x;d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %compute bottom out_disp_y            
            y_disp=[in_fic_bottom_y;d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_y=a*dis_x^2+b*dis_x+c;  
            %compute y direction derivative
            out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
        elseif(out_boun_2==2)
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %left situation
            j=in_fic_index(global_x-1,global_y);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_left_x=in_fic_for_out_ux(j,1);
                in_fic_left_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==2)
                        in_fic_left_x=in_fic_for_out_ux(j,k);
                        in_fic_left_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute left out_disp_x
            x_disp=[in_fic_left_x;d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_x=a*dis_y^2+b*dis_y+c;
            %compute left out_disp_y
            y_disp=[in_fic_left_y;d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_y=a*dis_y^2+b*dis_y+c;
            %right situation
            j=in_fic_index(global_x+1,global_y);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_right_x=in_fic_for_out_ux(j,1);
                in_fic_right_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==2)
                        in_fic_right_x=in_fic_for_out_ux(j,k);
                        in_fic_right_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute right out_disp_x
            x_disp=[in_fic_right_x;d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %compute right out_disp_y
            y_disp=[in_fic_right_y;d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,1)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,1)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
        end
    elseif(out_boun_1==2)
        %left situation:('L')
        
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper situation
            %ux
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            %upper_2 situation
            %ux
            x_disp=[d_d(global_x-1,global_y+2,1);d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+2,2);d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
            out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
        %bottom situatio:('B')
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right situation
            %ux
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            %right2 situation
            %ux
            x_disp=[d_d(global_x+2,global_y-1,1);d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+2,global_y-1,2);d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
            %}
            %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper situation
            %ux
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            %bottom situation
            %ux
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,2);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
        %bottom situatio:('B')
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right situation
            %ux
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            %left situation
            %ux
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,1);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,2)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,2)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
            %}
    elseif(out_boun_1==4)
        if(out_boun_2==2)
            %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper points, real points
            %compute upper out_disp_x
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %compute upper out_disp_y
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            % bottom situation
            j=in_fic_index(global_x,global_y-1);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_bottom_x=in_fic_for_out_ux(j,1);
                in_fic_bottom_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==1)
                        in_fic_bottom_x=in_fic_for_out_ux(j,k);
                        in_fic_bottom_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute bottom out_disp_x
            x_disp=[in_fic_bottom_x;d_d(global_x,global_y-1,1);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %compute bottom out_disp_y            
            y_disp=[in_fic_bottom_y;d_d(global_x,global_y-1,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_y=a*dis_x^2+b*dis_x+c;  
            %compute y direction derivative
            out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y); 
            %}
            
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %compute upper out_disp_y
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            % upper2 situation
            x_disp=[d_d(global_x-1,global_y+2,1);d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
            %compute upper out_disp_y
            y_disp=[d_d(global_x-1,global_y+2,2);d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
            %compute y direction derivative
            out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
            %}
        elseif(out_boun_2==1)
            %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %left situation
            j=in_fic_index(global_x-1,global_y);
            if(in_fic_for_out_ux(j,2)==0)
                in_fic_left_x=in_fic_for_out_ux(j,1);
                in_fic_left_y=in_fic_for_out_uy(j,1);
            elseif(in_fic_for_out_ux(j,2)~=0)
                for k=1:1:2
                    if(in_dir_index(j,k)==2)
                        in_fic_left_x=in_fic_for_out_ux(j,k);
                        in_fic_left_y=in_fic_for_out_uy(j,k);
                    end
                end
            end
            %compute left out_disp_x
            x_disp=[in_fic_left_x;d_d(global_x-1,global_y,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_x=a*dis_y^2+b*dis_y+c;
            %compute left out_disp_y
            y_disp=[in_fic_left_y;d_d(global_x-1,global_y,2);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_y=a*dis_y^2+b*dis_y+c;
            %right situation,real point
            %compute right out_disp_x
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %compute right out_disp_y
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,1)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,1)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
            %}
            
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right2 situation
            x_disp=[d_d(global_x+2,global_y-1,1);d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_x=a*dis_y^2+b*dis_y+c;
            %compute right out_disp_y
            y_disp=[d_d(global_x+2,global_y-1,2);d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_y=a*dis_y^2+b*dis_y+c;
            %right situation,real point
            %compute right out_disp_x
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %compute right out_disp_y
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
            %}
        end
    elseif(out_boun_1==5)
        if(out_boun_2==2)
        %left situation:('L')
        
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            j_1=in_fic_index(global_x,global_y+1);
            j_2=in_fic_index(global_x,global_y+2);
            if(point_material(global_x-1,global_y+2)==-1)
                    %ux
                x_disp=[in_fic_for_out_ux(j_1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                y_disp=[in_fic_for_out_uy(j_1,1);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                if(in_fic_for_out_ux(j_2,2)==0)
                    out_upper2_x=in_fic_for_out_ux(j_2,1);
                    out_upper2_y=in_fic_for_out_uy(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==1)
                            out_upper2_x=in_fic_for_out_ux(j_2,k);
                            out_upper2_y=in_fic_for_out_uy(j_2,k);
                        end
                    end
                end
                    %ux
                x_disp=[out_upper2_x;d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                y_disp=[out_upper2_y;d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
                out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
                out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
            elseif(point_material(global_x-1,global_y+2)==1)
                for k=1:1:2
                    if(in_dir_index(j_1,k)==1)
                        out_upper_x=in_fic_for_out_ux(j_1,k);
                        out_upper_y=in_fic_for_out_uy(j_1,k);
                    end
                end
                    %ux
                x_disp=[out_upper_x;d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                y_disp=[out_upper_y;d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                    %ux
                x_disp=[d_d(global_x-1,global_y+2,1);d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                y_disp=[d_d(global_x-1,global_y+2,2);d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
                syms a b c;
                eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
                out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
                out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
            end
        %}
        %bottom situatio:('B')
        %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,2);d_d(global_x+1,global_y-1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,2);d_d(global_x+1,global_y-1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            j=in_fic_index(global_x,global_y+1);
            x_disp=[in_fic_for_out_ux(j,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            y_disp=[in_fic_for_out_uy(j,1);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
            %}
            %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right situation
            %ux
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
        %left situation
            %ux
            j=in_fic_index(global_x,global_y+1);
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,1);in_fic_for_out_ux(j,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,1);in_fic_for_out_uy(j,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,2)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,2)=(out_right_disp_y-out_left_disp_y)/(2*d_x); 
            %}
            
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right situation
            %ux
            x_disp=[d_d(global_x+1,global_y-1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+1,global_y-1,2);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            %right2 situation
            %ux
            x_disp=[d_d(global_x+2,global_y-1,1);d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_x=a*dis_y^2+b*dis_y+c;
            %uy
            y_disp=[d_d(global_x+2,global_y-1,2);d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right2_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);   
          %}
        elseif(out_boun_2==1)
        %left situation:('L')
            %{
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper situation
            %ux
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            %bottom situation
            %ux
            j=in_fic_index(global_x+1,global_y);
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,2);in_fic_for_out_ux(j,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,2);in_fic_for_out_uy(j,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_bottom_disp_y=a*dis_x^2+b*dis_x+c;
            out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
            %}
            
            
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %upper situation
            %ux
            x_disp=[d_d(global_x-1,global_y+1,1);d_d(global_x,global_y+1,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+1,2);d_d(global_x,global_y+1,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper_disp_y=a*dis_x^2+b*dis_x+c;
            %upper_2 situation
            %ux
            x_disp=[d_d(global_x-1,global_y+2,1);d_d(global_x,global_y+2,1);d_d(global_x+1,global_y+2,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
            %uy
            y_disp=[d_d(global_x-1,global_y+2,2);d_d(global_x,global_y+2,2);d_d(global_x+1,global_y+2,2)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
            out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
            out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
            %}
        %bottom situatio:('B')
        %{
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,1);d_d(global_x-1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,1);d_d(global_x-1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_left_disp_y=a*dis_y^2+b*dis_y+c;
            j=in_fic_index(global_x+1,global_y);
            x_disp=[in_fic_for_out_ux(j,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_x=a*dis_y^2+b*dis_y+c;
            y_disp=[in_fic_for_out_uy(j,1);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            out_right_disp_y=a*dis_y^2+b*dis_y+c;
            out_deri_x_disp_x(i,2)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
            out_deri_x_disp_y(i,2)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
            %}
        
        
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %right situation
            j_1=in_fic_index(global_x+1,global_y);
            j_2=in_fic_index(global_x+2,global_y);
            if(point_material(global_x+2,global_y-1)==-1)
                x_disp=[in_fic_for_out_ux(j_1,1);d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[in_fic_for_out_uy(j_1,1);d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_y=a*dis_y^2+b*dis_y+c;
                if(in_fic_for_out_ux(j_2,2)==0)
                    in_right2_x=in_fic_for_out_ux(j_2,1);
                    in_right2_y=in_fic_for_out_uy(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==2)
                            in_right2_x=in_fic_for_out_ux(j_2,k);
                            in_right2_y=in_fic_for_out_uy(j_2,k);
                        end
                    end
                end
                x_disp=[in_right2_x;d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right2_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[in_right2_y;d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right2_disp_y=a*dis_y^2+b*dis_y+c;
                out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
                out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
            elseif(point_material(global_x+2,global_y-1)==1)
                for k=1:1:2
                    if(in_dir_index(j_1,k)==2)
                        in_right_x=in_fic_for_out_ux(j_1,k);
                        in_right_y=in_fic_for_out_uy(j_1,k);
                    end
                end
                x_disp=[in_right_x;d_d(global_x+1,global_y,1);d_d(global_x+1,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[in_right_y;d_d(global_x+1,global_y,2);d_d(global_x+1,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_y=a*dis_y^2+b*dis_y+c;
                x_disp=[d_d(global_x+2,global_y-1,1);d_d(global_x+2,global_y,1);d_d(global_x+2,global_y+1,1)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right2_disp_x=a*dis_y^2+b*dis_y+c;
                y_disp=[d_d(global_x+2,global_y-1,2);d_d(global_x+2,global_y,2);d_d(global_x+2,global_y+1,2)];
                syms a b c;
                eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right2_disp_y=a*dis_y^2+b*dis_y+c;
                out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
                out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
            end 
            %}
        end
    end
end