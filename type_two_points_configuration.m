function[]=type_two_points_configuration(size_type_two,d_x,d_y,r)
    global disp_gradient_x
    global disp_gradient_y
    global in_deri_x_disp_x
    global in_deri_x_disp_y
    global in_deri_y_disp_x
    global in_deri_y_disp_y
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global in_cross_boun_mark
    global out_fic_index
    global type_two_points
    global d_d
    
    for k=1:1:size_type_two
        i=type_two_points(k,1);
        j=type_two_points(k,2);
        k_t=type_two_points(k,3);
        disp_gradient_y(i+1,j,1)=(out_fic_for_in_ux(k_t,1)-d_d(i+1,j-1,1))/(2*d_y);
        disp_gradient_y(i+1,j,2)=(out_fic_for_in_uy(k_t,1)-d_d(i+1,j-1,2))/(2*d_y);
        disp_gradient_x(i,j+1,1)=(out_fic_for_in_ux(k_t,1)-d_d(i-1,j+1,1))/(2*d_x);
        disp_gradient_x(i,j+1,2)=(out_fic_for_in_uy(k_t,1)-d_d(i-1,j+1,2))/(2*d_x);
        %upper point
        t_1=out_fic_index(i,j+1);
        boun_1=in_cross_boun_mark(i,j+1,1);
        if(boun_1==5)
            global_x=i;
            global_y=j+1;
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %ux direction
            x_disp=[d_d(global_x-1,global_y,1);d_d(global_x,global_y,1);out_fic_for_in_ux(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_x_disp_x(t_1,1)=2*a*dis_x+b;
            %uy direction
            y_disp=[d_d(global_x-1,global_y,2);d_d(global_x,global_y,2);out_fic_for_in_uy(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_x_disp_y(t_1,1)=2*a*dis_x+b;
        elseif((boun_1==3)||(boun_1==4))
            global_x=i;
            global_y=j+1;
            dis_y=(global_y-1)*d_y;
            dis_x=sqrt(r^2-dis_y^2);
            %ux direction
            x_disp=[d_d(global_x-1,global_y,1);d_d(global_x,global_y,1);out_fic_for_in_ux(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==x_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==x_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_x_disp_x(t_1,1)=2*a*dis_x+b;
            %uy direction
            y_disp=[d_d(global_x-1,global_y,2);d_d(global_x,global_y,2);out_fic_for_in_uy(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_x-2)*d_x)^2+b*((global_x-2)*d_x)+c==y_disp(1);
            eqn(2)=a*((global_x-1)*d_x)^2+b*((global_x-1)*d_x)+c==y_disp(2);
            eqn(3)=a*((global_x)*d_x)^2+b*((global_x)*d_x)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_x_disp_y(t_1,1)=2*a*dis_x+b; 
        end
        %right point
        t_2=out_fic_index(i+1,j);
        boun_1=in_cross_boun_mark(i+1,j,1);
        if(boun_1==5)
            global_x=i+1;
            global_y=j;
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %ux direction
            x_disp=[d_d(global_x,global_y-1,1);d_d(global_x,global_y,1);out_fic_for_in_ux(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_y_disp_x(t_2,1)=2*a*dis_y+b;
            %uy direction
            y_disp=[d_d(global_x,global_y-1,2);d_d(global_x,global_y,2);out_fic_for_in_uy(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_y_disp_y(t_2,1)=2*a*dis_y+b;
        elseif((boun_1==3)||(boun_1==4))
            global_x=i+1;
            global_y=j;
            dis_x=(global_x-1)*d_x;
            dis_y=sqrt(r^2-dis_x^2);
            %ux direction
            x_disp=[d_d(global_x,global_y-1,1);d_d(global_x,global_y,1);out_fic_for_in_ux(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==x_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==x_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==x_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_y_disp_x(t_2,2)=2*a*dis_y+b;
            %uy direction
            y_disp=[d_d(global_x,global_y-1,2);d_d(global_x,global_y,2);out_fic_for_in_uy(k_t,1)];
            syms a b c;
            eqn(1)=a*((global_y-2)*d_y)^2+b*((global_y-2)*d_y)+c==y_disp(1);
            eqn(2)=a*((global_y-1)*d_y)^2+b*((global_y-1)*d_y)+c==y_disp(2);
            eqn(3)=a*((global_y)*d_y)^2+b*((global_y)*d_y)+c==y_disp(3);
            [a,b,c]=solve(eqn,[a,b,c]);
            in_deri_y_disp_y(t_2,2)=2*a*dis_y+b;
        end
    end
end