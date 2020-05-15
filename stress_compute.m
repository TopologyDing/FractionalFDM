function[]=stress_compute(r,d_x,d_y,high_x,high_y,max_x,max_y,miu_aluminum,lambda_aluminum,miu_steel,lambda_steel,size_out_fic_for_in,size_in_fic_for_out)
    global stress_x
    global stress_y
    global stress_xy
    global displacement
    global sec_order_frac_deri_x_disp_x
    global sec_order_frac_deri_x_disp_y
    global sec_order_frac_deri_y_disp_x
    global sec_order_frac_deri_y_disp_y
    global boun_frac_deri_x_disp_x
    global boun_frac_deri_x_disp_y
    global boun_frac_deri_y_disp_x
    global boun_frac_deri_y_disp_y
    global out_deri_x_disp_x
    global out_deri_x_disp_y
    global out_deri_y_disp_x
    global out_deri_y_disp_y
    global disp_gradient_x
    global disp_gradient_y
    global point_material
    global in_dir_index
    global out_dir_index
    global in_index_x
    global in_index_y
    global out_index_x
    global out_index_y
    global x_position
    global y_position
    deri_x_disp_x=zeros(max_x,max_y);
    deri_x_disp_y=zeros(max_x,max_y);
    deri_y_disp_x=zeros(max_x,max_y);
    deri_y_disp_y=zeros(max_x,max_y);
    
    for i=1:1:high_x
        for j=1:1:high_y
            if(sec_order_frac_deri_x_disp_x(i,j)~=0)
				symbolic_variables=symvar(sec_order_frac_deri_x_disp_x(i,j));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(sec_order_frac_deri_x_disp_x(i,j),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_x_disp_x(i,j)=deri_x_disp_x(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(sec_order_frac_deri_x_disp_y(i,j));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(sec_order_frac_deri_x_disp_y(i,j),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_x_disp_y(i,j)=deri_x_disp_y(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(sec_order_frac_deri_y_disp_x(i,j));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(sec_order_frac_deri_y_disp_x(i,j),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_y_disp_x(i,j)=deri_y_disp_x(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(sec_order_frac_deri_y_disp_y(i,j));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(sec_order_frac_deri_y_disp_y(i,j),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_y_disp_y(i,j)=deri_y_disp_y(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
            end
        end
    end
    for i=1:1:max_x
        for j=1:1:max_y
            if(point_material(i,j)>0)
				symbolic_variables=symvar(disp_gradient_x(i,j,1));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(disp_gradient_x(i,j,1),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_x_disp_x(i,j)=deri_x_disp_x(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(disp_gradient_x(i,j,2));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(disp_gradient_x(i,j,2),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_x_disp_y(i,j)=deri_x_disp_y(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(disp_gradient_y(i,j,1));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(disp_gradient_y(i,j,1),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_y_disp_x(i,j)=deri_y_disp_x(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
				symbolic_variables=symvar(disp_gradient_y(i,j,2));
				number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(disp_gradient_y(i,j,2),symbolic_variables));
				%analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    deri_y_disp_y(i,j)=deri_y_disp_y(i,j)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
            end
        end
    end
    
    
    for i=1:1:max_x
        for j=1:1:max_y
            if(point_material(i,j)<0)
                stress_x(i,j)=(2*miu_aluminum+lambda_aluminum)*(deri_x_disp_x(i,j))+lambda_aluminum*(deri_y_disp_y(i,j));
                stress_y(i,j)=(2*miu_aluminum+lambda_aluminum)*(deri_y_disp_y(i,j))+lambda_aluminum*(deri_x_disp_x(i,j));
                stress_xy(i,j)=miu_aluminum*(deri_x_disp_y(i,j)+deri_y_disp_x(i,j));
                x_position(i,j)=i;
                y_position(i,j)=j;
            else
                stress_x(i,j)=(2*miu_steel+lambda_steel)*(deri_x_disp_x(i,j))+lambda_steel*(deri_y_disp_y(i,j));
                stress_y(i,j)=(2*miu_steel+lambda_steel)*(deri_y_disp_y(i,j))+lambda_steel*(deri_x_disp_x(i,j));
                stress_xy(i,j)=miu_steel*(deri_x_disp_y(i,j)+deri_y_disp_x(i,j));
                x_position(i,j)=i;
                y_position(i,j)=j;
            end
        end
    end
    
    inboun_deri_x_disp_x=zeros(size_out_fic_for_in,2);
    inboun_deri_x_disp_y=zeros(size_out_fic_for_in,2);
    inboun_deri_y_disp_x=zeros(size_out_fic_for_in,2);
    inboun_deri_y_disp_y=zeros(size_out_fic_for_in,2);
    for i=1:1:size_out_fic_for_in
        for k=1:1:2
            fic_direction=out_dir_index(i,k);
            if(~fic_direction)
                continue;
            else
                symbolic_variables=symvar(boun_frac_deri_x_disp_x(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(boun_frac_deri_x_disp_x(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    inboun_deri_x_disp_x(i,k)=inboun_deri_x_disp_x(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(boun_frac_deri_x_disp_y(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(boun_frac_deri_x_disp_y(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    inboun_deri_x_disp_y(i,k)=inboun_deri_x_disp_y(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(boun_frac_deri_y_disp_x(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(boun_frac_deri_y_disp_x(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    inboun_deri_y_disp_x(i,k)=inboun_deri_y_disp_x(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(boun_frac_deri_y_disp_y(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(boun_frac_deri_y_disp_y(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    inboun_deri_y_disp_y(i,k)=inboun_deri_y_disp_y(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
            end
        end
    end
    outboun_deri_x_disp_x=zeros(size_in_fic_for_out,2);
    outboun_deri_x_disp_y=zeros(size_in_fic_for_out,2);
    outboun_deri_y_disp_x=zeros(size_in_fic_for_out,2);
    outboun_deri_y_disp_y=zeros(size_in_fic_for_out,2);
    for i=1:1:size_in_fic_for_out
        for k=1:1:2
            fic_direction=in_dir_index(i,k);
            if(~fic_direction)
                continue;
            else
                symbolic_variables=symvar(out_deri_x_disp_x(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(out_deri_x_disp_x(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    outboun_deri_x_disp_x(i,k)=outboun_deri_x_disp_x(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(out_deri_x_disp_y(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(out_deri_x_disp_y(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    outboun_deri_x_disp_y(i,k)=outboun_deri_x_disp_y(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(out_deri_y_disp_x(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(out_deri_y_disp_x(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    outboun_deri_y_disp_x(i,k)=outboun_deri_y_disp_x(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
                symbolic_variables=symvar(out_deri_y_disp_y(i,k));
                number_of_variable=length(symbolic_variables);
                coeff_variable=flip(coeffs(out_deri_y_disp_y(i,k),symbolic_variables));
                %analyse variables in each equilibrium expression
                for p=1:1:number_of_variable
                    str_variable=char(symbolic_variables(p));
                    length_of_string=length(str_variable);
                    %analyse each string
                    number_of_underline=0;
                    place=zeros(2,1);
                    for q=1:1:length_of_string
                        if(str_variable(q)=='_')
                            number_of_underline=number_of_underline+1;
                            place(number_of_underline)=q;
                            if(number_of_underline==2)
                                break;
                            else
                                continue;
                            end
                        end
                    end
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
                        sub_num_3=1;
                    else
                        sub_num_3=2;
                    end
                    string_1=str_variable(place(1)+1:1:place(2)-1);
                    string_2=str_variable(place(2)+1:1:length_of_string);
                    sub_num_1=str2double(string_1);
                    sub_num_2=str2double(string_2); 
                    outboun_deri_y_disp_y(i,k)=outboun_deri_y_disp_y(i,k)+coeff_variable(p)*displacement(sub_num_1,sub_num_2,sub_num_3);
                end
            end
        end
    end
    
    layer_size_1=ceil(size_out_fic_for_in/max_x);
    layer_size_2=ceil(size_in_fic_for_out/max_x);
    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        if(i~=max_x)
            rem_1=mod(i,max_x);
        else
            rem_1=max_x;
        end
        layer_1=fix(i/(max_x))+1;
        for k=1:1:2
            fic_direction=out_dir_index(i,k);
            if(fic_direction==1)
                stress_x(rem_1,max_y+layer_1)=(2*miu_aluminum+lambda_aluminum)*(inboun_deri_x_disp_x(i,k))+lambda_aluminum*(inboun_deri_y_disp_y(i,k));
                stress_y(rem_1,max_y+layer_1)=(2*miu_aluminum+lambda_aluminum)*(inboun_deri_y_disp_y(i,k))+lambda_aluminum*(inboun_deri_x_disp_x(i,k));
                stress_xy(rem_1,max_y+layer_1)=miu_aluminum*(inboun_deri_x_disp_y(i,k)+inboun_deri_y_disp_x(i,k));
                x_posi=sqrt(r^2-((global_y-1)*d_y)^2);
                x_position(rem_1,max_y+layer_1)=x_posi/d_x+1;
                y_position(rem_1,max_y+layer_1)=global_y;         
            elseif(fic_direction==2)
                stress_x(rem_1,max_y+layer_1+layer_size_1)=(2*miu_aluminum+lambda_aluminum)*(inboun_deri_x_disp_x(i,k))+lambda_aluminum*(inboun_deri_y_disp_y(i,k));
                stress_y(rem_1,max_y+layer_1+layer_size_1)=(2*miu_aluminum+lambda_aluminum)*(inboun_deri_y_disp_y(i,k))+lambda_aluminum*(inboun_deri_x_disp_x(i,k));
                stress_xy(rem_1,max_y+layer_1+layer_size_1)=miu_aluminum*(inboun_deri_x_disp_y(i,k)+inboun_deri_y_disp_x(i,k));
                y_posi=sqrt(r^2-((global_x-1)*d_x)^2);
                x_position(rem_1,max_y+layer_1+layer_size_1)=global_x;
                y_position(rem_1,max_y+layer_1+layer_size_1)=y_posi/d_y+1;
            end
        end
    end
    for i=1:1:size_in_fic_for_out
        global_x=in_index_x(i);
        global_y=in_index_y(i);
        if(i~=max_x)
            rem_2=mod(i,max_x);
        else
            rem_2=max_x;
        end
        layer_2=fix(i/(max_x))+1;
        for k=1:1:2
            fic_direction=in_dir_index(i,k);
            if(fic_direction==1)
                stress_x(rem_2,max_y+2*layer_size_1+layer_2)=(2*miu_steel+lambda_steel)*(outboun_deri_x_disp_x(i,k))+lambda_steel*(outboun_deri_y_disp_y(i,k));
                stress_y(rem_2,max_y+2*layer_size_1+layer_2)=(2*miu_steel+lambda_steel)*(outboun_deri_y_disp_y(i,k))+lambda_steel*(outboun_deri_x_disp_x(i,k));
                stress_xy(rem_2,max_y+2*layer_size_1+layer_2)=miu_steel*(outboun_deri_x_disp_y(i,k)+outboun_deri_y_disp_x(i,k));
                x_posi=sqrt(r^2-((global_y-1)*d_y)^2);
                x_position(rem_2,max_y+2*layer_size_1+layer_2)=x_posi/d_x+1;
                y_position(rem_2,max_y+2*layer_size_1+layer_2)=global_y;
            elseif(fic_direction==2)
                stress_x(rem_2,max_y+2*layer_size_1+layer_size_2+layer_2)=(2*miu_steel+lambda_steel)*(outboun_deri_x_disp_x(i,k))+lambda_steel*(outboun_deri_y_disp_y(i,k));
                stress_y(rem_2,max_y+2*layer_size_1+layer_size_2+layer_2)=(2*miu_steel+lambda_steel)*(outboun_deri_y_disp_y(i,k))+lambda_steel*(outboun_deri_x_disp_x(i,k));
                stress_xy(rem_2,max_y+2*layer_size_1+layer_size_2+layer_2)=miu_steel*(outboun_deri_x_disp_y(i,k)+outboun_deri_y_disp_x(i,k));
                y_posi=sqrt(r^2-((global_x-1)*d_x)^2);
                x_position(rem_2,max_y+2*layer_size_1+layer_size_2+layer_2)=global_x;
                y_position(rem_2,max_y+2*layer_size_1+layer_size_2+layer_2)=y_posi/d_y+1;
            end
        end
    end
    
    x_position(x_position==0)=NaN;
    y_position(y_position==0)=NaN;
    stress_x(stress_x==0)=NaN;
    stress_y(stress_y==0)=NaN;
    stress_xy(stress_xy==0)=NaN;
    
    
    
end