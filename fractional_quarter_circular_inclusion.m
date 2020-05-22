function[]=fractional_quarter_circular_inclusion()
    dbstop if error

    %material properties
    v_steel=0.28883542806974205;
    miu_steel=8.237068454836098e10;
    lambda_steel=1.1266838804655762e11;
    v_steel=0.3;
    miu_steel=2;
    lambda_steel=2*miu_steel*v_steel/(1-2*v_steel);
    %v_aluminum=0.28883542806974205;
    %miu_aluminum=8.237068454836098e10;
    %lambda_aluminum=1.1266838804655762e11;
    v_aluminum=0.33124956931271565;
    miu_aluminum=2.596732215483388e10;
    lambda_aluminum=5.097269526934259e10;
    v_aluminum=0.4;
    miu_aluminum=1;
    lambda_aluminum=2*miu_aluminum*v_aluminum/(1-2*v_aluminum);
%     v_steel=0.33124956931271565;
%     lambda_steel=5.097269526934259e10;
%     miu_steel=2.596732215483388e10;
    %geometry configuration
    global r
    l_x=1.5;
    l_y=1.5;
    r=0.875;
    center_x=0;
    center_y=0;
    
    %nonlocal configuration
    global alpha
    alpha=1;
    GM_alpha=gamma(2-alpha);
    nonlocalrange=0;
    left_x=nonlocalrange;
    right_x=nonlocalrange;
    down_y=nonlocalrange;
    upper_y=nonlocalrange;
    nonlocal_size(1)=left_x;
    nonlocal_size(2)=right_x;
    nonlocal_size(3)=down_y;
    nonlocal_size(4)=upper_y;
    
    %discretization configuration
    global d_x
    global d_y
    d_x=0.06;
    d_y=0.06;
    max_x=l_x/d_x+1;
    max_y=l_y/d_y+1;
    
    global d_d
    d_d=sym(zeros(max_x,max_y,2));
    Ux=sym('Ux_',[max_x max_y]);
    Uy=sym('Uy_',[max_x max_y]);
    for i=1:1:max_x
        for j=1:1:max_y
            d_d(i,j,1)=Ux(i,j);
            d_d(i,j,2)=Uy(i,j);
        end
    end
    
    global point_material
    point_material=zeros(max_x,max_y);
    for i=1:1:max_x
        for j=1:1:max_y
            if(((i-1)*d_x)^2+((j-1)*d_x)^2<r^2)
                point_material(i,j)=-2;
            else
                point_material(i,j)=2;
            end
        end
    end
    high_x=ceil(r/d_x)+1;
    high_y=ceil(r/d_y)+1;
    
    boun_mark=zeros(high_x+1,high_y+1);
    for i=1:1:high_x+1
        for j=1:1:high_y+1
            real_x=(i-1)*d_x-center_x;
            real_y=(j-1)*d_y-center_y;
            if(((real_x)^2+(real_y)^2)<r^2)
                if((((real_x+d_x)^2+(real_y)^2)>r^2)||(((real_x)^2+(real_y+d_y)^2)>r^2)...
                        ||(((real_x-d_x)^2+(real_y)^2)>r^2)||(((real_x)^2+(real_y-d_y)^2)>r^2)...
                        ||(((real_x+d_x)^2+(real_y+d_y)^2)>r^2)||(((real_x+d_x)^2+(real_y-d_y)^2)>r^2)...
                        ||(((real_x-d_x)^2+(real_y+d_y)^2)>r^2)||(((real_x-d_x)^2+(real_y-d_y)^2)>r^2))
                    boun_mark(i,j)=-1;
                    point_material(i,j)=-1;
                end
            else
                if((((real_x+d_x)^2+(real_y)^2)<r^2)||(((real_x)^2+(real_y+d_y)^2)<r^2)...
                        ||(((real_x-d_x)^2+(real_y)^2)<r^2)||(((real_x)^2+(real_y-d_y)^2)<r^2)...
                        ||(((real_x+d_x)^2+(real_y+d_y)^2)<r^2)||(((real_x+d_x)^2+(real_y-d_y)^2)<r^2)...
                        ||(((real_x-d_x)^2+(real_y+d_y)^2)<r^2)||(((real_x-d_x)^2+(real_y-d_y)^2)<r^2))
                    boun_mark(i,j)=1;
                    point_material(i,j)=1;
                end
            end
        end
    end
    
    % in and out boundary configuration
    global in_cross_boun_mark
    global out_cross_boun_mark
    in_cross_boun_mark=zeros(high_x,high_y,2);
	out_cross_boun_mark=zeros(high_x,high_y,2);
    for i=2:1:high_x
        for j=2:1:high_y
            if(boun_mark(i,j)==1)
                if((boun_mark(i-1,j)==-1)&&...
                    (boun_mark(i-1,j+1)==-1)&&...
                    (boun_mark(i-1,j-1)==-1))
                    out_cross_boun_mark(i,j,1)=1;
                    out_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j-1)==-1)&&...
                    (boun_mark(i,j-1)==-1)&&...
                    (boun_mark(i+1,j-1)==-1))
                    out_cross_boun_mark(i,j,1)=1;
                    out_cross_boun_mark(i,j,2)=2;
                end
                if((boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i,j-1)==-1)&&...
                    (boun_mark(i-1,j)==-1))
                    out_cross_boun_mark(i,j,1)=2;
                    out_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i,j-1)==1)&&...
                    (boun_mark(i-1,j)==1))
                    out_cross_boun_mark(i,j,1)=3;
                    out_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j)==1)&&...
                    (boun_mark(i-1,j)==1)&&...
                    (boun_mark(i,j-1)==-1)&&...
                    (boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i+1,j-1)==1))
                    out_cross_boun_mark(i,j,1)=4;
                    out_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i,j+1)==1)&&...
                    (boun_mark(i,j-1)==1)&&...
                    (boun_mark(i-1,j)==-1)&&...
                    (boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i-1,j+1)==1))
                    out_cross_boun_mark(i,j,1)=4;
                    out_cross_boun_mark(i,j,2)=2;
                end
                if((boun_mark(i+1,j)==1)&&...
                    (boun_mark(i-1,j)==-1)&&...
                    (boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i,j-1)==-1)&&...
                    (boun_mark(i+1,j-1)==-1))
                    out_cross_boun_mark(i,j,1)=5;
                    out_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i,j+1)==1)&&...
                    (boun_mark(i-1,j+1)==-1)&&...
                    (boun_mark(i-1,j)==-1)&&...
                    (boun_mark(i-1,j-1)==-1)&&...
                    (boun_mark(i,j-1)==-1))
                    out_cross_boun_mark(i,j,1)=5;
                    out_cross_boun_mark(i,j,2)=2;
                end
            elseif(boun_mark(i,j)==-1)
                if((boun_mark(i+1,j)==1)&&...
                    (boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i+1,j-1)==1))
                        in_cross_boun_mark(i,j,1)=1;
                        in_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i,j+1)==1)&&...
                    (boun_mark(i+1,j+1)==1))
                        in_cross_boun_mark(i,j,1)=1;
                        in_cross_boun_mark(i,j,2)=2;
                end
                if((boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i,j+1)==-1)&&...
                    (boun_mark(i+1,j)==-1))
                        in_cross_boun_mark(i,j,1)=2;
                        in_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i,j+1)==1)&&...
                    (boun_mark(i+1,j)==1))
                        in_cross_boun_mark(i,j,1)=3;
                        in_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j)==1)&&...
                    (boun_mark(i-1,j+1)==1)&&...
                    (boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i,j+1)==1)&&...
                    (boun_mark(i-1,j)==-1))
                        in_cross_boun_mark(i,j,1)=4;
                        in_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i,j+1)==1)&&...
                    (boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i+1,j)==1)&&...
                    (boun_mark(i+1,j-1)==1)&&...
                    (boun_mark(i,j-1)==-1))
                        in_cross_boun_mark(i,j,1)=4;
                        in_cross_boun_mark(i,j,2)=2;
                end
                if((boun_mark(i,j+1)==1)&&...
                    (boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i-1,j+1)==-1)&&...
                    (boun_mark(i-1,j)==-1)&&...
                    (boun_mark(i+1,j)==-1))
                        in_cross_boun_mark(i,j,1)=5;
                        in_cross_boun_mark(i,j,2)=1;
                end
                if((boun_mark(i+1,j)==1)&&...
                    (boun_mark(i+1,j+1)==1)&&...
                    (boun_mark(i,j+1)==-1)&&...
                    (boun_mark(i,j-1)==-1)&&...
                    (boun_mark(i+1,j-1)==-1))
                        in_cross_boun_mark(i,j,1)=5;
                        in_cross_boun_mark(i,j,2)=2;
                end
            end		
        end
    end
    local_point=point_material(1:1:high_x,1:1:high_y);
	%{
    for i=1:1:high_x
        for j=1:1:high_y
            if(local_point(i,j)==1)
                scatter3(i,j,local_point(i,j),'filled','red');
                hold on;
            elseif(local_point(i,j)==-1)
                scatter3(i,j,local_point(i,j),'filled','blue');
                hold on;
            else
                scatter3(i,j,local_point(i,j),'filled','green');
                hold on;
            end
        end
    end
	%}
	in_fic_for_out=zeros(high_x,high_y);
	out_fic_for_in=zeros(high_x,high_y);
    
    for i=2:1:high_x
        for j=2:1:high_y
            if(out_cross_boun_mark(i,j,1)==1)
                if(out_cross_boun_mark(i,j,2)==1)
                    in_fic_for_out(i,j)=1;
                elseif(out_cross_boun_mark(i,j,2)==2)
                    in_fic_for_out(i,j)=2;
                end
            elseif(out_cross_boun_mark(i,j,1)==2)
                in_fic_for_out(i,j)=3;
            elseif(out_cross_boun_mark(i,j,1)==3)
                in_fic_for_out(i,j)=4;
            elseif(out_cross_boun_mark(i,j,1)==4)
                if(out_cross_boun_mark(i,j,2)==1)
                    in_fic_for_out(i,j)=2;
                elseif(out_cross_boun_mark(i,j,2)==2)
                    in_fic_for_out(i,j)=1;
                end
            elseif(out_cross_boun_mark(i,j,1)==5)
                in_fic_for_out(i,j)=3;
            end
        end
    end
    
    for i=2:1:high_x
        for j=2:1:high_y
            if(in_cross_boun_mark(i,j,1)==1)
                if(in_cross_boun_mark(i,j,2)==1)
                    out_fic_for_in(i,j)=1;
                elseif(in_cross_boun_mark(i,j,2)==2)
                    out_fic_for_in(i,j)=2;
                end
            elseif(in_cross_boun_mark(i,j,1)==2)
                out_fic_for_in(i,j)=3;
            elseif(in_cross_boun_mark(i,j,1)==3)
                out_fic_for_in(i,j)=4;
            elseif(in_cross_boun_mark(i,j,1)==4)
                out_fic_for_in(i,j)=4;
            elseif(in_cross_boun_mark(i,j,1)==5)
                if(in_cross_boun_mark(i,j,2)==1)
                    out_fic_for_in(i,j)=2;
                elseif(in_cross_boun_mark(i,j,2)==2)
                    out_fic_for_in(i,j)=1;
                end
            end
        end
    end
    %implement first point and last point into boundary mark
    % in_cross_boun_mark(i,j,1)=-1: boundary point on left or down side
    % in_cross_boun_mark(i,j,2)=1: on left boundary
    % in_cross_boun_mark(i,j,2)=2: on right boundary
    in_cross_boun_mark(1,high_y-1,1)=-1;
    in_cross_boun_mark(1,high_y-1,2)=2;
    in_cross_boun_mark(high_x-1,1,1)=-1;
    in_cross_boun_mark(high_x-1,1,2)=1;
    out_cross_boun_mark(1,high_y,1)=-1;
    out_cross_boun_mark(1,high_y,2)=2;
    out_cross_boun_mark(high_x,1,1)=-1;
    out_cross_boun_mark(high_x,1,2)=1;
	
    %first_deri_matrix=zeros(max_x,max_y,2,3,3);
    % in and out fic point configuration
    global in_fic_for_out_ux
    global in_fic_for_out_uy
    global out_fic_for_in_ux
    global out_fic_for_in_uy
    global in_dir_index
    global out_dir_index
    global in_index_x
    global in_index_y
    global out_index_x
    global out_index_y
	in_fic_for_out_ux=sym(zeros(high_x*5,3));
	in_fic_for_out_uy=sym(zeros(high_x*5,3));
	in_dir_index=zeros(high_x*5,3);
	in_index_x=zeros(high_x*5);
	in_index_y=zeros(high_x*5);
	index=1;
	out_fic_for_in_ux=sym(zeros(high_x*5,3));
	out_fic_for_in_uy=sym(zeros(high_x*5,3));
	out_dir_index=zeros(high_x*5,3);
	out_index_x=zeros(high_x*5);
	out_index_y=zeros(high_x*5);
	index_1=1;
    
    in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(1)+'_'+string(high_y)+'_B');
    in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(1)+'_'+string(high_y)+'_B');
    in_dir_index(index,1)=2;
    in_index_x(index)=1;
    in_index_y(index)=high_y;
    index=index+1;
    out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(1)+'_'+string(high_y-1)+'_T');
    out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(1)+'_'+string(high_y-1)+'_T');
    out_dir_index(index_1,1)=2;
    out_index_x(index_1)=1;
    out_index_y(index_1)=high_y-1;
    index_1=index_1+1;
    for i=2:1:high_x
        for j=high_y:-1:2
            if(in_fic_for_out(i,j)==1)
                in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_L');
                in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_L');
                in_dir_index(index,1)=1;
                in_index_x(index)=i;
                in_index_y(index)=j;
                index=index+1;
            elseif(in_fic_for_out(i,j)==2)
                in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_B');
                in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_B');
                in_dir_index(index,1)=2;
                in_index_x(index)=i;
                in_index_y(index)=j;
                index=index+1;
            elseif(in_fic_for_out(i,j)==3)
                in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_L');
                in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_L');
                in_dir_index(index,1)=1;
                in_fic_for_out_ux(index,2)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_B');
                in_fic_for_out_uy(index,2)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_B');
                in_dir_index(index,2)=2;
                in_fic_for_out_ux(index,3)=3;
                in_fic_for_out_uy(index,3)=3;
                in_index_x(index)=i;
                in_index_y(index)=j;
                index=index+1;
            elseif(in_fic_for_out(i,j)==4)
                in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_LB');
                in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_LB');
                in_index_x(index)=i;
                in_index_y(index)=j;
                index=index+1;
            end
            if(out_fic_for_in(i,j)==1)
                out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_R');
                out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_R');
                out_dir_index(index_1,1)=1;
                out_index_x(index_1)=i;
                out_index_y(index_1)=j;
                index_1=index_1+1;
            elseif(out_fic_for_in(i,j)==2)
                out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_T');
                out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_T');
                out_dir_index(index_1,1)=2;
                out_index_x(index_1)=i;
                out_index_y(index_1)=j;
                index_1=index_1+1;
            elseif(out_fic_for_in(i,j)==3)
                out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_RT');
                out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_RT');
                out_index_x(index_1)=i;
                out_index_y(index_1)=j;
                index_1=index_1+1;
            elseif(out_fic_for_in(i,j)==4)
                out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_R');
                out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_R');
                out_dir_index(index_1,1)=1;
                out_fic_for_in_ux(index_1,2)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_T');
                out_fic_for_in_uy(index_1,2)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_T');
                out_dir_index(index_1,2)=2;
                out_fic_for_in_ux(index_1,3)=4;
                out_fic_for_in_uy(index_1,3)=4;
                out_index_x(index_1)=i;
                out_index_y(index_1)=j;
                index_1=index_1+1;
            end
        end
    end
    
    in_fic_for_out_ux(index,1)=sym('in_fic_ux_'+string(high_x)+'_'+string(1)+'_L');
    in_fic_for_out_uy(index,1)=sym('in_fic_uy_'+string(high_x)+'_'+string(1)+'_L');
    in_dir_index(index,1)=1;
    in_index_x(index)=high_x;
    in_index_y(index)=1;
    out_fic_for_in_ux(index_1,1)=sym('out_fic_ux_'+string(high_x-1)+'_'+string(1)+'_R');
    out_fic_for_in_uy(index_1,1)=sym('out_fic_uy_'+string(high_x-1)+'_'+string(1)+'_R');
    out_dir_index(index_1,1)=1;
    out_index_x(index_1)=high_x-1;
    out_index_y(index_1)=1;
    
    for i=1:1:length(in_fic_for_out_ux)
        if(in_fic_for_out_ux(i,1)==0)
            in_fic_for_out_ux=in_fic_for_out_ux(1:1:i-1,:);
            in_fic_for_out_uy=in_fic_for_out_uy(1:1:i-1,:);
            in_index_x=in_index_x(1:1:i-1);
            in_index_y=in_index_y(1:1:i-1);
            in_dir_index=in_dir_index(1:1:i-1,:);
            break;
        end
    end
    for i=1:1:length(out_fic_for_in_ux)
        if(out_fic_for_in_ux(i,1)==0)
            out_fic_for_in_ux=out_fic_for_in_ux(1:1:i-1,:);
            out_fic_for_in_uy=out_fic_for_in_uy(1:1:i-1,:);
            out_index_x=out_index_x(1:1:i-1);
            out_index_y=out_index_y(1:1:i-1);
            out_dir_index=out_dir_index(1:1:i-1,:);
            break;
        end
    end
    
    global in_fic_index
    global out_fic_index
	in_fic_index=zeros(high_x+1,high_y+1);
	out_fic_index=zeros(high_x+1,high_y+1);
    for k=1:1:length(in_fic_for_out_ux)
        for i=1:1:high_x
            for j=1:1:high_y
                if((in_index_x(k)==i)&&(in_index_y(k)==j))
                    in_fic_index(i,j)=k;
                end
            end
        end
    end
	
    for k=1:1:length(out_fic_for_in_ux)
        for i=1:1:high_x
            for j=1:1:high_y
                if((out_index_x(k)==i)&&(out_index_y(k)==j))
                    out_fic_index(i,j)=k;
                end
            end
        end
    end
    
    for i=1:length(in_fic_for_out_ux)
        if(in_fic_for_out_ux(i,3)==3)
            g_x=in_index_x(i);
            g_y=in_index_y(i);
            if(point_material(g_x+1,g_y-1)==1)
                R_D_x=d_d(g_x+1,g_y-1,1);
                R_D_y=d_d(g_x+1,g_y-1,2);
            else
                index_R_D=in_fic_index(g_x+1,g_y);
                for k=1:2
                    if(in_dir_index(index_R_D,k)==2)
                        R_D_x=in_fic_for_out_ux(index_R_D,k);
                        R_D_y=in_fic_for_out_uy(index_R_D,k);
                        break;
                    end
                end
            end
            if(point_material(g_x-1,g_y+1)==1)
                L_U_x=d_d(g_x-1,g_y+1,1);
                L_U_y=d_d(g_x-1,g_y+1,2);
            else
                index_L_U=in_fic_index(g_x,g_y+1);
                for k=1:2
                    if(in_dir_index(index_L_U,k)==1)
                        L_U_x=in_fic_for_out_ux(index_L_U,k);
                        L_U_y=in_fic_for_out_uy(index_L_U,k);
                        break;
                    end
                end
            end
            in_fic_for_out_ux(i,3)=(4*in_fic_for_out_ux(i,1)+4*in_fic_for_out_ux(i,2)-4*d_d(g_x,g_y,1)+...
                (d_d(g_x+1,g_y+1,1)-L_U_x-R_D_x))/3;
            in_fic_for_out_uy(i,3)=(4*in_fic_for_out_uy(i,1)+4*in_fic_for_out_uy(i,2)-4*d_d(g_x,g_y,2)+...
                (d_d(g_x+1,g_y+1,2)-L_U_y-R_D_y))/3;
        end
    end
    for i=1:length(out_fic_for_in_ux)
        if(out_fic_for_in_ux(i,3)==4)
            g_x=out_index_x(i);
            g_y=out_index_y(i);
            if(point_material(g_x+1,g_y-1)==-1)
                R_D_x=d_d(g_x+1,g_y-1,1);
                R_D_y=d_d(g_x+1,g_y-1,2);
            else
                index_R_D=out_fic_index(g_x,g_y-1);
                for k=1:2
                    if(out_dir_index(index_R_D,k)==1)
                        R_D_x=out_fic_for_in_ux(index_R_D,k);
                        R_D_y=out_fic_for_in_uy(index_R_D,k);
                        break;
                    end
                end
            end
            if(point_material(g_x-1,g_y+1)==-1)
                L_U_x=d_d(g_x-1,g_y+1,1);
                L_U_y=d_d(g_x-1,g_y+1,2);
            else
                index_L_U=out_fic_index(g_x-1,g_y);
                for k=1:2
                    if(out_dir_index(index_L_U,k)==2)
                        L_U_x=out_fic_for_in_ux(index_L_U,k);
                        L_U_y=out_fic_for_in_uy(index_L_U,k);
                        break;
                    end
                end
            end
            out_fic_for_in_ux(i,3)=(4*out_fic_for_in_ux(i,1)+4*out_fic_for_in_ux(i,2)-4*d_d(g_x,g_y,1)+...
                (-L_U_x-R_D_x+d_d(g_x-1,g_y-1,1)))/3;
            out_fic_for_in_uy(i,3)=(4*out_fic_for_in_uy(i,1)+4*out_fic_for_in_uy(i,2)-4*d_d(g_x,g_y,2)+...
                (-L_U_y-R_D_y+d_d(g_x-1,g_y-1,2)))/3;
        end
    end    
    %}
    
    %% 1st order displacement gradient configuration
    global disp_gradient_x
    global disp_gradient_y
    disp_gradient_x=sym(zeros(max_x,max_y,2));
    disp_gradient_y=sym(zeros(max_x,max_y,2));
    for i=1:1:max_x
        for j=1:1:max_y
            t=point_material(i,j);
            if(t==1)
                boun_point_type=in_fic_for_out(i,j);
                results=integer_displacement_gradient(i,j,t,boun_point_type,max_x,max_y);
            elseif(t==-1)
                boun_point_type=out_fic_for_in(i,j);
                results=integer_displacement_gradient(i,j,t,boun_point_type,max_x,max_y);
            else
                results=integer_displacement_gradient(i,j,t,0,max_x,max_y);
            end
            disp_gradient_x(i,j,1)=results(1);
            disp_gradient_x(i,j,2)=results(2);
            disp_gradient_y(i,j,1)=results(3);
            disp_gradient_y(i,j,2)=results(4);
        end
    end
    
    global one_order_disp_gradient_x
    global one_order_disp_gradient_y
    one_order_disp_gradient_x=sym(zeros(high_x,high_y,2));
    one_order_disp_gradient_y=sym(zeros(high_x,high_y,2));
    one_order_precision_disp_gradient(high_x,high_y,d_x,d_y);
    %% 2nd order interpolation for computing boundary displacement and gradient
    global out_sur_nor_dir_x
    global out_sur_nor_dir_y
    global in_deri_x_disp_x
    global in_deri_x_disp_y
    global in_deri_y_disp_x
    global in_deri_y_disp_y
    global out_deri_x_disp_x
    global out_deri_x_disp_y
    global out_deri_y_disp_x
    global out_deri_y_disp_y
    global in_disp_x
    global in_disp_y
    global out_disp_x
    global out_disp_y
    size_in_fic_for_out=length(in_fic_for_out_ux);
    size_out_fic_for_in=length(out_fic_for_in_ux);  
    
    out_sur_nor_dir_x=zeros(size_out_fic_for_in,3);
    out_sur_nor_dir_y=zeros(size_out_fic_for_in,3);
    
    in_disp_x=sym(zeros(size_out_fic_for_in,3));
    in_disp_y=sym(zeros(size_out_fic_for_in,3));
    out_disp_x=sym(zeros(size_in_fic_for_out,3));
    out_disp_y=sym(zeros(size_in_fic_for_out,3));

    in_deri_x_disp_x=sym(zeros(size_out_fic_for_in,3));
    in_deri_x_disp_y=sym(zeros(size_out_fic_for_in,3));
    in_deri_y_disp_x=sym(zeros(size_out_fic_for_in,3));
    in_deri_y_disp_y=sym(zeros(size_out_fic_for_in,3));

    out_deri_x_disp_x=sym(zeros(size_in_fic_for_out,3));
    out_deri_x_disp_y=sym(zeros(size_in_fic_for_out,3));
    out_deri_y_disp_x=sym(zeros(size_in_fic_for_out,3));
    out_deri_y_disp_y=sym(zeros(size_in_fic_for_out,3));
    
    for i=1:1:size_in_fic_for_out
        global_x=in_index_x(i);
        global_y=in_index_y(i);
        out_boun_1=out_cross_boun_mark(global_x,global_y,1);
        out_boun_2=out_cross_boun_mark(global_x,global_y,2);
        out_boun_disp_grad(i,out_boun_1,out_boun_2,global_x,global_y,size_in_fic_for_out);
    end
    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        in_boun_1=in_cross_boun_mark(global_x,global_y,1);
        in_boun_2=in_cross_boun_mark(global_x,global_y,2);
        in_boun_disp_grad(i,in_boun_1,in_boun_2,global_x,global_y,size_out_fic_for_in);
    end
    in_disp_x=vpa(in_disp_x);
    in_disp_y=vpa(in_disp_y);
    in_deri_x_disp_x=vpa(in_deri_x_disp_x);
    in_deri_x_disp_y=vpa(in_deri_x_disp_y);
    in_deri_y_disp_x=vpa(in_deri_y_disp_x);
    in_deri_y_disp_y=vpa(in_deri_y_disp_y);
    out_disp_x=vpa(out_disp_x);
    out_disp_y=vpa(out_disp_y);
    out_deri_x_disp_x=vpa(out_deri_x_disp_x);
    out_deri_x_disp_y=vpa(out_deri_x_disp_y);
    out_deri_y_disp_x=vpa(out_deri_y_disp_x);
    out_deri_y_disp_y=vpa(out_deri_y_disp_y);
    %% nonlocal range configuration for boundary
    %point position for x direction derivative:
    %point_posi_x_dir(1):x coord
    %point_posi_x_dir(2):y coord
    %point position for y direction derivative:
    %point_posi_y_dir(1):x coord
    %point_posi_y_dir(2):y coord
    point_posi_x_dir=ones(high_x,high_y,2)*-1;
    point_posi_y_dir=ones(high_x,high_y,2)*-1;
    %boundary points position
    global boun_x
    global boun_y
    boun_x=zeros(high_y,2);
    boun_y=zeros(high_x,2);
    global boun_points
    boun_points=zeros(high_x,high_y);
    left_boun_x_posi=zeros(high_y,2);
    right_boun_x_posi=zeros(high_y,2);
    down_boun_y_posi=zeros(high_x,2);
    upper_boun_y_posi=zeros(high_x,2);
    
    for i=1:1:high_x
        distance_x=(i-1)*d_x;
        boun_y(i,1)=0;
        if(distance_x>r)
            boun_y(i,1)=-1;
            boun_y(i,2)=-1;
        else
            distance_y=sqrt(r^2-distance_x^2);
            boun_y(i,2)=distance_y;
        end
    end
    for j=1:1:high_y
        boun_x(j,1)=0;
        distance_y=(j-1)*d_y;
        if(distance_y>r)
            boun_x(j,1)=-1;
            boun_x(j,2)=-1;
        else
            distance_x=sqrt(r^2-distance_y^2);
            boun_x(j,2)=distance_x;
        end
    end
    for i=1:1:high_x
        for j=1:1:high_y
            if(point_material(i,j)==-2)
                point_posi_x_dir(i,j,1)=(i-1)*d_x;
                point_posi_x_dir(i,j,2)=(j-1)*d_y;
                point_posi_y_dir(i,j,1)=(i-1)*d_x;
                point_posi_y_dir(i,j,2)=(j-1)*d_y;
            elseif(point_material(i,j)==-1)
                point_posi_x_dir(i,j,1)=(i-1)*d_x;
                point_posi_x_dir(i,j,2)=(j-1)*d_y;
                point_posi_y_dir(i,j,1)=(i-1)*d_x;
                point_posi_y_dir(i,j,2)=(j-1)*d_y; 
                if(point_material(i+1,j)==1)
                    boun_points(i+1,j)=1;
                    point_posi_x_dir(i+1,j,1)=boun_x(j,2);
                    point_posi_x_dir(i+1,j,2)=(j-1)*d_y;
                    right_boun_x_posi(j,1)=i+1;
                    right_boun_x_posi(j,2)=j;
                end
                if(point_material(i,j+1)==1)
                    point_posi_y_dir(i,j+1,1)=(i-1)*d_x;
                    point_posi_y_dir(i,j+1,2)=boun_y(i,2);
                    upper_boun_y_posi(i,1)=i;
                    upper_boun_y_posi(i,2)=j+1;
                end
            end
            left_boun_x_posi(j,1)=1;
            left_boun_x_posi(j,2)=j;
            down_boun_y_posi(i,1)=i;
            down_boun_y_posi(i,2)=1;
        end
    end
    
    
    %boundary info configuration
    % boun_info:
    % boun_info(:,:,1):position of nonlocal boundary on each direction
    % boun_info(:,:,2):length of nonlocal boundary on each direction
    % boun_info(:,:,3):if nonlocal boundary on the geometric boundary
    % point_posi:
    % point_posi(:,:,1):x position of the point
    % point_posi(:,:,2):y position of the point
    % point_posi(:,:,3):point on x direction or y direction
    % point_posi(:,:,3)==1:point on x direction
    % point_posi(:,:,3)==2:point on y direction
    global boun_left_info
    global boun_right_info
    global boun_down_info
    global boun_upper_info
    global boun_point_posi
    boun_left_info=ones(size_out_fic_for_in,2,3)*-1;
    boun_right_info=ones(size_out_fic_for_in,2,3)*-1;
    boun_down_info=ones(size_out_fic_for_in,2,3)*-1;
    boun_upper_info=ones(size_out_fic_for_in,2,3)*-1;
    boun_point_posi=zeros(size_out_fic_for_in,2,3);
    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        in_boun(1)=in_cross_boun_mark(global_x,global_y,1);
        in_boun(2)=in_cross_boun_mark(global_x,global_y,2);
        boun_info(in_boun,i,global_x,global_y,nonlocal_size);
    end
   
    %% nonlocal range configuration for points not on the boundary
    global left_info
    global right_info
    global down_info
    global upper_info
    left_info=ones(high_x,high_y,3)*-1;
    right_info=ones(high_x,high_y,3)*-1;
    down_info=ones(high_x,high_y,3)*-1;
    upper_info=ones(high_x,high_y,3)*-1;
    for i=1:1:high_x
        for j=1:1:high_y
            % x direction configuration
            left_point_x=1;
            right_point_x=right_boun_x_posi(j,1);
            if((i>left_point_x+left_x)&&(i<right_point_x+1))
                left_info(i,j,2)=left_x;
                left_info(i,j,1)=(i-left_x-1)*d_x;
            elseif((i>left_point_x-1)&&(i<right_point_x+1))
                left_info(i,j,3)=1;
                left_info(i,j,2)=i-left_point_x;
                left_info(i,j,1)=boun_x(j,1);
            end
            if((i<right_point_x-right_x)&&(i>left_point_x-1))
                right_info(i,j,2)=right_x;
                right_info(i,j,1)=(i+right_x-1)*d_x;
            elseif((i<right_point_x+1)&&(i>left_point_x-1))
                right_info(i,j,3)=1;
                right_info(i,j,2)=right_point_x-i;
                right_info(i,j,1)=boun_x(j,2);
            end 
            %y direction configuration
            down_point_y=1;
            upper_point_y=upper_boun_y_posi(i,2);
            if((j>down_point_y+down_y)&&(j<upper_point_y+1))
                down_info(i,j,2)=down_y;
                down_info(i,j,1)=(j-down_y-1)*d_y;
            elseif((j>down_point_y-1)&&(j<upper_point_y+1))
                down_info(i,j,3)=1;
                down_info(i,j,2)=j-down_point_y;
                down_info(i,j,1)=boun_y(i,1);
            end
            if((j<upper_point_y-upper_y)&&(j>down_point_y-1))
                upper_info(i,j,2)=upper_y;
                upper_info(i,j,1)=(j+upper_y-1)*d_y;
            elseif((j<upper_point_y+1)&&(j>down_point_y-1))
                upper_info(i,j,3)=1;
                upper_info(i,j,2)=upper_point_y-j;
                upper_info(i,j,1)=boun_y(i,2);
            end
        end
    end
    
    %% nonlocal derivative computation for boundary points
    global boun_frac_deri_x_disp_x
    global boun_frac_deri_x_disp_y
    global boun_frac_deri_y_disp_x
    global boun_frac_deri_y_disp_y
    global boun_frac_deri_x_disp_x_l
    global boun_frac_deri_x_disp_x_r
    global boun_frac_deri_x_disp_y_l
    global boun_frac_deri_x_disp_y_r
    global boun_frac_deri_y_disp_x_d
    global boun_frac_deri_y_disp_x_u
    global boun_frac_deri_y_disp_y_d
    global boun_frac_deri_y_disp_y_u
    boun_frac_deri_x_disp_x=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_x_disp_y=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_x=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_y=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_x_disp_x_l=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_x_disp_y_l=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_x_d=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_y_d=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_x_disp_x_r=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_x_disp_y_r=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_x_u=sym(zeros(size_out_fic_for_in,2));
    boun_frac_deri_y_disp_y_u=sym(zeros(size_out_fic_for_in,2));

    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        for k=1:1:2
            if(out_dir_index(i,k)~=0)
                result_l=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_left_info(i,k,:),i,k,global_x,global_y,1);
                result_r=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_right_info(i,k,:),i,k,global_x,global_y,2);
                result_d=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_down_info(i,k,:),i,k,global_x,global_y,3);
                result_u=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_upper_info(i,k,:),i,k,global_x,global_y,4);
                if(boun_left_info(i,k,2)==0)
                    boun_frac_deri_x_disp_x_l(i,k)=result_l(1)/GM_alpha;
                    boun_frac_deri_x_disp_y_l(i,k)=result_l(2)/GM_alpha;
                else
                    boun_frac_deri_x_disp_x_l(i,k)=result_l(1)*(boun_point_posi(i,k,1)-boun_left_info(i,k,1))^(alpha-1);
                    boun_frac_deri_x_disp_y_l(i,k)=result_l(2)*(boun_point_posi(i,k,1)-boun_left_info(i,k,1))^(alpha-1);
                end
                if(boun_right_info(i,k,2)==0)
                    boun_frac_deri_x_disp_x_r(i,k)=result_r(1)/GM_alpha;
                    boun_frac_deri_x_disp_y_r(i,k)=result_r(2)/GM_alpha;
                else
                    boun_frac_deri_x_disp_x_r(i,k)=result_r(1)*(boun_right_info(i,k,1)-boun_point_posi(i,k,1))^(alpha-1);
                    boun_frac_deri_x_disp_y_r(i,k)=result_r(2)*(boun_right_info(i,k,1)-boun_point_posi(i,k,1))^(alpha-1);
                end
                if(boun_down_info(i,k,2)==0)
                    boun_frac_deri_y_disp_x_d(i,k)=result_d(1)/GM_alpha;
                    boun_frac_deri_y_disp_y_d(i,k)=result_d(2)/GM_alpha;
                else
                    boun_frac_deri_y_disp_x_d(i,k)=result_d(1)*(boun_point_posi(i,k,2)-boun_down_info(i,k,1))^(alpha-1);
                    boun_frac_deri_y_disp_y_d(i,k)=result_d(2)*(boun_point_posi(i,k,2)-boun_down_info(i,k,1))^(alpha-1);
                end
                if(boun_upper_info(i,k,2)==0)
                    boun_frac_deri_y_disp_x_u(i,k)=result_u(1)/GM_alpha;
                    boun_frac_deri_y_disp_y_u(i,k)=result_u(2)/GM_alpha;
                else
                    boun_frac_deri_y_disp_x_u(i,k)=result_u(1)*(boun_upper_info(i,k,1)-boun_point_posi(i,k,2))^(alpha-1);
                    boun_frac_deri_y_disp_y_u(i,k)=result_u(2)*(boun_upper_info(i,k,1)-boun_point_posi(i,k,2))^(alpha-1);
                end
            else
                continue;
            end
            boun_frac_deri_x_disp_x(i,k)=vpa((GM_alpha*(boun_frac_deri_x_disp_x_l(i,k)+boun_frac_deri_x_disp_x_r(i,k)))/2);
            boun_frac_deri_x_disp_y(i,k)=vpa((GM_alpha*(boun_frac_deri_x_disp_y_l(i,k)+boun_frac_deri_x_disp_y_r(i,k)))/2);
            boun_frac_deri_y_disp_x(i,k)=vpa((GM_alpha*(boun_frac_deri_y_disp_x_d(i,k)+boun_frac_deri_y_disp_x_u(i,k)))/2);
            boun_frac_deri_y_disp_y(i,k)=vpa((GM_alpha*(boun_frac_deri_y_disp_y_d(i,k)+boun_frac_deri_y_disp_y_u(i,k)))/2);
        end
    end
    %}
    global sym_boun_point_posi
    global sym_left_info
    global sym_right_info
    global sym_down_info
    global sym_upper_info
    global finite_dis
    sym_boun_point_posi=zeros(size_out_fic_for_in,2,3);
    sym_left_info=ones(size_out_fic_for_in,2,3)*-1;
    sym_right_info=ones(size_out_fic_for_in,2,3)*-1;
    sym_down_info=ones(size_out_fic_for_in,2,3)*-1;
    sym_upper_info=ones(size_out_fic_for_in,2,3)*-1;
    finite_dis=zeros(size_out_fic_for_in,2);
    compute_boun_corre_point_info(size_out_fic_for_in,d_x,d_y,nonlocal_size,r);
    %% nonlocal derivative computation for inner points
    %{
    for i=1:1:high_x
        for j=1:1:high_y
            if(point_material(i,j)<0)
                x_posi=(i-1)*d_x;
                y_posi=(j-1)*d_y;
                results=Caputo_derivative(left_info(i,j,:),i,j,1);
                lleft_x=results(1);
                lleft_y=results(2);
                results=Caputo_derivative(right_info(i,j,:),i,j,2);
                lright_x=results(1);
                lright_y=results(2);
                results=Caputo_derivative(down_info(i,j,:),i,j,3);
                ldown_x=results(1);
                ldown_y=results(2);
                results=Caputo_derivative(upper_info(i,j,:),i,j,4);
                lupper_x=results(1);
                lupper_y=results(2);
                if(left_info(i,j,2)==0)
                    frac_deri_x_disp_x_l(i,j)=lleft_x/gamma(2-alpha);
                    frac_deri_x_disp_y_l(i,j)=lleft_y/gamma(2-alpha);
                else
                    frac_deri_x_disp_x_l(i,j)=lleft_x*(x_posi-left_info(i,j,1))^(alpha-1);
                    frac_deri_x_disp_y_l(i,j)=lleft_y*(x_posi-left_info(i,j,1))^(alpha-1);
                end
                if(right_info(i,j,2)==0)
                    frac_deri_x_disp_x_r(i,j)=lright_x/gamma(2-alpha);
                    frac_deri_x_disp_y_r(i,j)=lright_y/gamma(2-alpha);
                else
                    frac_deri_x_disp_x_r(i,j)=lright_x*(right_info(i,j,1)-x_posi)^(alpha-1);
                    frac_deri_x_disp_y_r(i,j)=lright_y*(right_info(i,j,1)-x_posi)^(alpha-1);
                end
                if(down_info(i,j,2)==0)
                    frac_deri_y_disp_x_d(i,j)=ldown_x/gamma(2-alpha);
                    frac_deri_y_disp_y_d(i,j)=ldown_y/gamma(2-alpha);
                else
                    frac_deri_y_disp_x_d(i,j)=ldown_x*(y_posi-down_info(i,j,1))^(alpha-1);
                    frac_deri_y_disp_y_d(i,j)=ldown_y*(y_posi-down_info(i,j,1))^(alpha-1);
                end
                if(upper_info(i,j,2)==0)
                    frac_deri_y_disp_x_u(i,j)=lupper_x/gamma(2-alpha);
                    frac_deri_y_disp_y_u(i,j)=lupper_y/gamma(2-alpha);
                else
                    frac_deri_y_disp_x_u(i,j)=lupper_x*(upper_info(i,j,1)-y_posi)^(alpha-1);
                    frac_deri_y_disp_y_u(i,j)=lupper_y*(upper_info(i,j,1)-y_posi)^(alpha-1);
                end
                frac_deri_x_disp_x(i,j)=vpa((gamma(2-alpha)*(frac_deri_x_disp_x_l(i,j)+frac_deri_x_disp_x_r(i,j)))/2);
                frac_deri_x_disp_y(i,j)=vpa((gamma(2-alpha)*(frac_deri_x_disp_y_l(i,j)+frac_deri_x_disp_y_r(i,j)))/2);
                frac_deri_y_disp_x(i,j)=vpa((gamma(2-alpha)*(frac_deri_y_disp_x_d(i,j)+frac_deri_y_disp_x_u(i,j)))/2);
                frac_deri_y_disp_y(i,j)=vpa((gamma(2-alpha)*(frac_deri_y_disp_y_d(i,j)+frac_deri_y_disp_y_u(i,j)))/2);
            end
        end
    end
    %}  
    %% boundary conditions
    
    
	in_fic_vari_sqn=zeros(2*length(out_deri_x_disp_x),1);
	in_fic_variables=sym(zeros(2*length(out_deri_x_disp_x),1));
	out_fic_variables=sym(zeros(2*length(out_deri_x_disp_x),1));
	out_fic_vari_sqn=zeros(2*length(out_deri_x_disp_x),1);
	index_in=1;
	index_in_sqn=1;
	index_out=1;
	index_out_sqn=1;
    for i=1:1:length(in_fic_for_out_ux)
        if(in_fic_for_out_ux(i,2)==0)
            str_variable=char(in_fic_for_out_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1)=='_')
                in_fic_variables(index_in)=in_fic_for_out_ux(i,1);
                in_fic_variables(index_in+1)=in_fic_for_out_uy(i,1);
                in_fic_vari_sqn(index_in_sqn)=i;
                index_in=index_in+2;
                index_in_sqn=index_in_sqn+1;
            else
                continue;
            end
        else
            in_fic_variables(index_in)=in_fic_for_out_ux(i,1);
            in_fic_variables(index_in+1)=in_fic_for_out_uy(i,1);
            in_fic_vari_sqn(index_in_sqn)=i;
            in_fic_variables(index_in+2)=in_fic_for_out_ux(i,2);
            in_fic_variables(index_in+3)=in_fic_for_out_uy(i,2);
            in_fic_vari_sqn(index_in_sqn+1)=i;
            index_in_sqn=index_in_sqn+2;
            index_in=index_in+4;
        end
    end
    for i=1:1:length(out_fic_for_in_ux)
        if(out_fic_for_in_ux(i,2)==0)
            str_variable=char(out_fic_for_in_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1)=='_')
                out_fic_variables(index_out)=out_fic_for_in_ux(i,1);
                out_fic_variables(index_out+1)=out_fic_for_in_uy(i,1);
                out_fic_vari_sqn(index_out_sqn)=i;
                index_out=index_out+2;
                index_out_sqn=index_out_sqn+1;
            else
                continue;
            end
        else
            out_fic_variables(index_out)=out_fic_for_in_ux(i,1);
            out_fic_variables(index_out+1)=out_fic_for_in_uy(i,1);
            out_fic_vari_sqn(index_out_sqn)=i;
            out_fic_variables(index_out+2)=out_fic_for_in_ux(i,2);
            out_fic_variables(index_out+3)=out_fic_for_in_uy(i,2);
            out_fic_vari_sqn(index_out_sqn+1)=i;
            index_out=index_out+4;
            index_out_sqn=index_out_sqn+2;
        end
    end
    for i=1:1:length(in_fic_variables)
        if(in_fic_variables(i)==0)
            in_fic_variables=in_fic_variables(1:1:i-1);
            break;
        end
    end
    for i=1:1:length(out_fic_variables)
        if(out_fic_variables(i)==0)
            out_fic_variables=out_fic_variables(1:1:i-1);
            break;
        end
    end
    for i=1:1:length(in_fic_vari_sqn)
        if(in_fic_vari_sqn(i)==0)
            in_fic_vari_sqn=in_fic_vari_sqn(1:1:i-1);
            break;
        end
    end
    for i=1:1:length(out_fic_vari_sqn)
        if(out_fic_vari_sqn(i)==0)
            out_fic_vari_sqn=out_fic_vari_sqn(1:1:i-1);
            break;
        end
    end
    
	out_sur_nor_dir_x=vpa(out_sur_nor_dir_x);
	out_sur_nor_dir_y=vpa(out_sur_nor_dir_y);
	out_sur_nor_dir_x=out_sur_nor_dir_x/r;
	out_sur_nor_dir_y=out_sur_nor_dir_y/r;
	M_steel=2*miu_steel+lambda_steel;
	M_aluminum=2*miu_aluminum+lambda_aluminum;
    %{
	syms sin_t cos_t;
    boun_traction_1=[M_aluminum*cos_t,-M_steel*cos_t,miu_aluminum*sin_t,-miu_steel*sin_t,...
        miu_aluminum*sin_t,-miu_steel*sin_t,M_aluminum*v_aluminum*cos_t/(1-v_aluminum),-lambda_steel*v_steel*cos_t/(1-v_steel)];
    boun_traction_2=[M_aluminum*v_aluminum*sin_t/(1-v_aluminum),-lambda_steel*v_steel*sin_t/(1-v_steel),miu_aluminum*cos_t,-miu_steel*cos_t,...
        miu_aluminum*cos_t,-miu_steel*cos_t,M_aluminum*sin_t,-M_steel*sin_t];
    boun_traction_1=[M_aluminum*cos_t,-M_steel*cos_t,miu_aluminum*sin_t,-miu_steel*sin_t,...
        miu_aluminum*sin_t,-miu_steel*sin_t,M_aluminum*v_aluminum*cos_t/(1-v_aluminum),-M_steel*v_steel*cos_t/(1-v_steel)];
    boun_traction_2=[M_aluminum*v_aluminum*sin_t/(1-v_aluminum),-M_steel*v_steel*sin_t/(1-v_steel),miu_aluminum*cos_t,-miu_steel*cos_t,...
        miu_aluminum*cos_t,-miu_steel*cos_t,M_aluminum*sin_t,-M_steel*sin_t];
    
    M_p=M_aluminum;
    M_m=M_steel;
    m_p=miu_aluminum;
    m_m=miu_steel;
    l_p=M_aluminum*v_aluminum/(1-v_aluminum);
    l_m=M_steel*v_steel/(1-v_steel);
    boun_1=[-M_p*cos_t^2-m_m*sin_t^2,M_m*cos_t^2+m_m*sin_t^2,(m_m-m_p)*sin_t*cos_t,0,...
        -(l_m+m_p)*sin_t*cos_t,(l_m+m_m)*cos_t*sin_t,(l_m-l_p)*cos_t^2,0];
    boun_2=[-(l_p+m_m)*sin_t*cos_t,(l_m+m_m)*cos_t*sin_t,(m_m-m_p)*cos_t^2,0,...
        -m_p*cos_t^2-M_m*sin_t^2,m_m*cos_t^2+M_m*sin_t^2,(M_m-M_p)*sin_t*cos_t,0];
    boun_3=[M_p*cos_t^2+m_p*sin_t^2,-(M_m*cos_t^2+m_p*sin_t^2),0,(m_p-m_m)*sin_t*cos_t,...
        (m_p+l_p)*sin_t*cos_t,-(m_m+l_p)*sin_t*cos_t,0,(l_p-l_m)*cos_t^2];
    boun_4=[(l_p+m_p)*sin_t*cos_t,-(l_m+m_p)*sin_t*cos_t,0,(m_p-m_m)*cos_t^2,...
        m_p*cos_t^2+M_p*sin_t^2,-m_m*cos_t^2-M_p*sin_t^2,0,(M_p-M_m)*sin_t*cos_t];
    boun_5=[(M_p-M_m)*sin_t*cos_t,0,m_p*sin_t^2+M_m*cos_t^2,-m_m*sin_t^2-M_m*cos_t^2,...
        (m_p-m_m)*sin_t^2,0,(l_p+m_m)*sin_t*cos_t,-(l_m+m_m)*sin_t*cos_t];
    boun_6=[(l_p-l_m)*sin_t^2,0,(m_p+l_m)*sin_t*cos_t,-(m_m+l_m)*sin_t*cos_t,...
        (m_p-m_m)*sin_t*cos_t,0,M_p*sin_t^2+m_m*cos_t^2,-M_m*sin_t^2-m_m*cos_t^2];
    boun_7=[0,(M_m-M_p)*sin_t*cos_t,-(m_p*sin_t^2+M_p*cos_t^2),m_m*sin_t^2+M_p*cos_t^2,...
        0,(m_m-m_p)*sin_t^2,-(l_p+m_p)*sin_t*cos_t,(l_m+m_p)*sin_t*cos_t];
    boun_8=[0,(l_m-l_p)*sin_t^2,-(m_p+l_p)*sin_t*cos_t,(m_m+l_p)*sin_t*cos_t,...
        0,(m_m-m_p)*sin_t*cos_t,-M_p*sin_t^2-m_p*cos_t^2,M_m*sin_t^2+m_p*cos_t^2];
    
	length_variables=length(out_fic_variables);
	disconti_boun_eqn=sym(zeros(length_variables*2,1));
	index_eqn=1;
	
    for i=1:1:length(out_fic_for_in_ux)
        x_coor=out_index_x(i);
        y_coor=out_index_y(i);
        if(out_fic_for_in_ux(i,2)==0)
            %str_variable=char(out_fic_for_in_ux(i,1));
            %length_str=length(str_variable);
            if(out_dir_index(i,1)==1)
                j=in_fic_index(x_coor+1,y_coor);
            elseif(out_dir_index(i,1)==2)
                j=in_fic_index(x_coor,y_coor+1);
            else
                continue;
            end
            num_boun_trac_1=subs(boun_traction_1,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            num_boun_trac_2=subs(boun_traction_2,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            if(in_fic_for_out_ux(j,2)==0)
                C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,1);...
                    boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,1)];
                disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,1)==0;
                disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,1)==0;
                disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                index_eqn=index_eqn+4;
            else
                if(out_dir_index(i,1)==1)
                    out_boun=out_cross_boun_mark(x_coor+1,y_coor,1);
                elseif(out_dir_index(i,1)==2)
                    out_boun=out_cross_boun_mark(x_coor,y_coor+1,1);
                else
                    continue;
                end
                for m=1:1:2
                    if(out_boun==5)
                        if(out_dir_index(i,1)==1)
                            if(in_dir_index(j,m)==1)
                            num_boun_1=subs(boun_1,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
                            num_boun_2=subs(boun_2,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
                            C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,1)==2)
                            if(in_dir_index(j,m)==2)
                            num_boun_1=subs(boun_5,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
                            num_boun_2=subs(boun_6,[sin_t,cos_t],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
                            C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        end
                    else
                        if(out_dir_index(i,1)==1)
                            if(in_dir_index(j,m)==1)
                            C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,1)==2)
                            if(in_dir_index(j,m)==2)
                            C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        end
                    end
                end
            end
        elseif(in_cross_boun_mark(x_coor,y_coor,1)==3)
            for k=1:1:2
                if(out_dir_index(i,k)==1)
                    j=in_fic_index(x_coor+1,y_coor);
                elseif(out_dir_index(i,k)==2)
                    j=in_fic_index(x_coor,y_coor+1);
                end
                num_boun_trac_1=subs(boun_traction_1,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                num_boun_trac_2=subs(boun_traction_2,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                if(in_fic_for_out_ux(j,2)==0)
                    C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,1);...
                        boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,1)];
                    disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,1)==0;
                    disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,1)==0;
                    disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                    disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                    index_eqn=index_eqn+4;
                else
                    for m=1:1:2
                        if(out_dir_index(i,k)==1)
                            if(in_dir_index(j,m)==1)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,k)==2)
                            if(in_dir_index(j,m)==2)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        end
                    end
                end
            end    
        elseif(in_cross_boun_mark(x_coor,y_coor,1)==4)
            for k=1:1:2
                if(out_dir_index(i,k)==1)
                    j=in_fic_index(x_coor+1,y_coor);
                elseif(out_dir_index(i,k)==2)
                    j=in_fic_index(x_coor,y_coor+1);
                end
                num_boun_trac_1=subs(boun_traction_1,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                num_boun_trac_2=subs(boun_traction_2,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                if(in_fic_for_out_ux(j,2)==0)
                    if(out_dir_index(i,k)==1)
                        num_boun_1=subs(boun_3,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                        num_boun_2=subs(boun_4,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                    C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,1);...
                        boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,1)];
                    disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,1)==0;
                    disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,1)==0;
                    disconti_boun_eqn(index_eqn+2)=num_boun_1*C==0;
                    disconti_boun_eqn(index_eqn+3)=num_boun_2*C==0;
                    index_eqn=index_eqn+4;
                    elseif(out_dir_index(i,k)==2)
                        num_boun_1=subs(boun_7,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                        num_boun_2=subs(boun_8,[sin_t,cos_t],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                    C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,1);...
                        boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,1)];
                    disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,1)==0;
                    disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,1)==0;
                    disconti_boun_eqn(index_eqn+2)=num_boun_1*C==0;
                    disconti_boun_eqn(index_eqn+3)=num_boun_2*C==0;
                    index_eqn=index_eqn+4;
                    end
                else
                    for m=1:1:2
                        if(out_dir_index(i,k)==1)
                            if(in_dir_index(j,m)==1)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,k)==2)
                            if(in_dir_index(j,m)==2)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        end
                    end
                end
            end 
        end
    end
    %}
    
	syms sin_theta cos_theta;
    %{
    boun_traction_1=[M_aluminum*cos_theta,-M_steel*cos_theta,miu_aluminum*sin_theta,-miu_steel*sin_theta,...
        miu_aluminum*sin_theta,-miu_steel*sin_theta,M_aluminum*v_aluminum*cos_theta/(1-v_aluminum),-M_steel*v_steel*cos_theta/(1-v_steel)];
    boun_traction_2=[M_aluminum*v_aluminum*sin_theta/(1-v_aluminum),-M_steel*v_steel*sin_theta/(1-v_steel),miu_aluminum*cos_theta,-miu_steel*cos_theta,...
        miu_aluminum*cos_theta,-miu_steel*cos_theta,M_aluminum*sin_theta,-M_steel*sin_theta];
    %}
    boun_traction_1=[M_aluminum*cos_theta,-M_steel*cos_theta,miu_aluminum*sin_theta,-miu_steel*sin_theta,...
        miu_aluminum*sin_theta,-miu_steel*sin_theta,lambda_aluminum*cos_theta,-lambda_steel*cos_theta];
    boun_traction_2=[lambda_aluminum*sin_theta,-lambda_steel*sin_theta,miu_aluminum*cos_theta,-miu_steel*cos_theta,...
        miu_aluminum*cos_theta,-miu_steel*cos_theta,M_aluminum*sin_theta,-M_steel*sin_theta];
    %{
    boun_traction_1=[M_aluminum*cos_theta,-M_steel*cos_theta,miu_aluminum*sin_theta,-miu_steel*sin_theta,...
        miu_aluminum*sin_theta,-miu_steel*sin_theta,M_aluminum*v_aluminum*cos_theta/(1-v_aluminum),-lambda_steel*v_steel*cos_theta/(1-v_steel)];
    boun_traction_2=[M_aluminum*v_aluminum*sin_theta/(1-v_aluminum),-lambda_steel*v_steel*sin_theta/(1-v_steel),miu_aluminum*cos_theta,-miu_steel*cos_theta,...
        miu_aluminum*cos_theta,-miu_steel*cos_theta,M_aluminum*sin_theta,-M_steel*sin_theta];
    %}
    
    length_variables=length(out_fic_variables);
	disconti_boun_eqn=sym(zeros(length_variables*2,1));
	index_eqn=1;
    
    for i=1:1:size_out_fic_for_in
        x_coor=out_index_x(i);
        y_coor=out_index_y(i);
        if(out_fic_for_in_ux(i,2)==0)
            if(out_dir_index(i,1)==1)
                j=in_fic_index(x_coor+1,y_coor);
            elseif(out_dir_index(i,1)==2)
                j=in_fic_index(x_coor,y_coor+1);
            else
                continue;
            end
            num_boun_trac_1=subs(boun_traction_1,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            num_boun_trac_2=subs(boun_traction_2,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            if(in_fic_for_out_ux(j,2)==0)
                C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,1);...
                    boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,1)];
                disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,1)==0;
                disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,1)==0;
                disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                index_eqn=index_eqn+4;
            else
                for m=1:1:2
                    if(out_dir_index(i,1)==1)
                        if(in_dir_index(j,m)==1)
                        C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                            boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                        disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                        disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                        disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                        disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                        index_eqn=index_eqn+4;
                        end
                    elseif(out_dir_index(i,1)==2)
                        if(in_dir_index(j,m)==2)
                        C=[boun_frac_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                            boun_frac_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                        disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                        disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                        disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                        disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                        index_eqn=index_eqn+4;
                        end
                    end
                end
            end
        else
            for k=1:1:2
                if(out_dir_index(i,k)==1)
                    j=in_fic_index(x_coor+1,y_coor);
                elseif(out_dir_index(i,k)==2)
                    j=in_fic_index(x_coor,y_coor+1);
                end
                num_boun_trac_1=subs(boun_traction_1,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                num_boun_trac_2=subs(boun_traction_2,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,k),out_sur_nor_dir_x(i,k)]);
                if(in_fic_for_out_ux(j,2)==0)
                    C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,1);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,1);...
                        boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,1);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,1)];
                    disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,1)==0;
                    disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,1)==0;
                    disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                    disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                    index_eqn=index_eqn+4;
                else
                    for m=1:1:2
                        if(out_dir_index(i,k)==1)
                            if(in_dir_index(j,m)==1)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,k)==2)
                            if(in_dir_index(j,m)==2)
                            C=[boun_frac_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);boun_frac_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                boun_frac_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);boun_frac_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        end
                    end
                end
            end
        end
    end
    %}
	for i=1:1:length(disconti_boun_eqn)
		if(disconti_boun_eqn(i)==0)
			disconti_boun_eqn=disconti_boun_eqn(1:1:i-1);
			break;
		end
	end
	disconti_boun_eqn=vpa(disconti_boun_eqn);
	%% solve disconti_boun_eqns
	fic_variables=[in_fic_variables;out_fic_variables];
	solutions=solve(disconti_boun_eqn,fic_variables);
	fields=fieldnames(solutions);
    
	for i=1:1:length(in_fic_vari_sqn)
		j=in_fic_vari_sqn(i);
		if(i==1)
			if(j~=in_fic_vari_sqn(i+1))
				in_fic_for_out_ux(j,1)=solutions.(fields{2*i-1});
				in_fic_for_out_uy(j,1)=solutions.(fields{2*i});
			elseif(j==in_fic_vari_sqn(i+1))
				in_fic_for_out_ux(j,1)=solutions.(fields{2*i-1});
				in_fic_for_out_uy(j,1)=solutions.(fields{2*i});
				in_fic_for_out_ux(j,2)=solutions.(fields{2*i+1});
				in_fic_for_out_uy(j,2)=solutions.(fields{2*i+2});
			end
		elseif(i<length(in_fic_vari_sqn))
			if(j==in_fic_vari_sqn(i-1))
				continue;
			else
				if(j~=in_fic_vari_sqn(i+1))
					in_fic_for_out_ux(j,1)=solutions.(fields{2*i-1});
					in_fic_for_out_uy(j,1)=solutions.(fields{2*i});
				elseif(j==in_fic_vari_sqn(i+1))
					in_fic_for_out_ux(j,1)=solutions.(fields{2*i-1});
					in_fic_for_out_uy(j,1)=solutions.(fields{2*i});
					in_fic_for_out_ux(j,2)=solutions.(fields{2*i+1});
					in_fic_for_out_uy(j,2)=solutions.(fields{2*i+2});
				end
			end
		else
			if(j==in_fic_vari_sqn(i-1))
				continue;
			else
				in_fic_for_out_ux(j,1)=solutions.(fields{2*i-1});
				in_fic_for_out_uy(j,1)=solutions.(fields{2*i});
			end
		end
	end

	for i=1:1:length(out_fic_vari_sqn)
		j=out_fic_vari_sqn(i);
		if(i==1)
			if(j~=out_fic_vari_sqn(i+1))
				out_fic_for_in_ux(j,1)=solutions.(fields{2*i+length_variables-1});
				out_fic_for_in_uy(j,1)=solutions.(fields{2*i+length_variables});
			elseif(j==out_fic_vari_sqn(i+1))
				out_fic_for_in_ux(j,1)=solutions.(fields{2*i+length_variables-1});
				out_fic_for_in_uy(j,1)=solutions.(fields{2*i+length_variables});
				out_fic_for_in_ux(j,2)=solutions.(fields{2*i+length_variables+1});
				out_fic_for_in_uy(j,2)=solutions.(fields{2*i+length_variables+2});
			end
		elseif(i<length(out_fic_vari_sqn))
			if(j==out_fic_vari_sqn(i-1))
				continue;
			else
				if(j~=out_fic_vari_sqn(i+1))
					out_fic_for_in_ux(j,1)=solutions.(fields{2*i+length_variables-1});
					out_fic_for_in_uy(j,1)=solutions.(fields{2*i+length_variables});
				elseif(j==out_fic_vari_sqn(i+1))
					out_fic_for_in_ux(j,1)=solutions.(fields{2*i+length_variables-1});
					out_fic_for_in_uy(j,1)=solutions.(fields{2*i+length_variables});
					out_fic_for_in_ux(j,2)=solutions.(fields{2*i+length_variables+1});
					out_fic_for_in_uy(j,2)=solutions.(fields{2*i+length_variables+2});
				end
			end
		else
			if(j==out_fic_vari_sqn(i-1))
				continue;
			else
				out_fic_for_in_ux(j,1)=solutions.(fields{2*i+length_variables-1});
				out_fic_for_in_uy(j,1)=solutions.(fields{2*i+length_variables});
			end
		end
	end
	
	for i=1:1:length(in_fic_for_out_ux)
		g_x=in_index_x(i);
		g_y=in_index_y(i);
		if(out_cross_boun_mark(g_x,g_y,1)==3)
            in_fic_for_out_ux(i,1)=(4*d_d(g_x-1,g_y,1)+4*d_d(g_x,g_y-1,1)-4*d_d(g_x,g_y,1)+...
                (d_d(g_x+1,g_y+1,1)-d_d(g_x+1,g_y-1,1)-d_d(g_x-1,g_y+1,1)))/3;
            in_fic_for_out_uy(i,1)=(4*d_d(g_x-1,g_y,2)+4*d_d(g_x,g_y-1,2)-4*d_d(g_x,g_y,2)+...
                (d_d(g_x+1,g_y+1,2)-d_d(g_x+1,g_y-1,2)-d_d(g_x-1,g_y+1,2)))/3;
		end
		if(in_fic_for_out_ux(i,3)~=0)
            g_x=in_index_x(i);
            g_y=in_index_y(i);
            if(point_material(g_x+1,g_y-1)==1)
                R_D_x=d_d(g_x+1,g_y-1,1);
                R_D_y=d_d(g_x+1,g_y-1,2);
            else
                index_R_D=in_fic_index(g_x+1,g_y);
                for k=1:2
                    if(in_dir_index(index_R_D,k)==2)
                        R_D_x=in_fic_for_out_ux(index_R_D,k);
                        R_D_y=in_fic_for_out_uy(index_R_D,k);
                        break;
                    end
                end
            end
            if(point_material(g_x-1,g_y+1)==1)
                L_U_x=d_d(g_x-1,g_y+1,1);
                L_U_y=d_d(g_x-1,g_y+1,2);
            else
                index_L_U=in_fic_index(g_x,g_y+1);
                for k=1:2
                    if(in_dir_index(index_L_U,k)==1)
                        L_U_x=in_fic_for_out_ux(index_L_U,k);
                        L_U_y=in_fic_for_out_uy(index_L_U,k);
                        break;
                    end
                end
            end
            in_fic_for_out_ux(i,3)=(4*in_fic_for_out_ux(i,1)+4*in_fic_for_out_ux(i,2)-4*d_d(g_x,g_y,1)+...
                (d_d(g_x+1,g_y+1,1)-L_U_x-R_D_x))/3;
            in_fic_for_out_uy(i,3)=(4*in_fic_for_out_uy(i,1)+4*in_fic_for_out_uy(i,2)-4*d_d(g_x,g_y,2)+...
                (d_d(g_x+1,g_y+1,2)-L_U_y-R_D_y))/3;
		end
	end
    for i=1:1:length(out_fic_for_in_ux)
        g_x=out_index_x(i);
        g_y=out_index_y(i);
        if(in_cross_boun_mark(g_x,g_y,1)==2)
            out_fic_for_in_ux(i,1)=(4*d_d(g_x+1,g_y,1)+4*d_d(g_x,g_y+1,1)-4*d_d(g_x,g_y,1)...
                -d_d(g_x-1,g_y+1,1)-d_d(g_x+1,g_y-1,1)+d_d(g_x-1,g_y-1,1))/3;
            out_fic_for_in_uy(i,1)=(4*d_d(g_x+1,g_y,2)+4*d_d(g_x,g_y+1,2)-4*d_d(g_x,g_y,2)...
                -d_d(g_x-1,g_y+1,2)-d_d(g_x+1,g_y-1,2)+d_d(g_x-1,g_y-1,2))/3;
        end
        if(out_fic_for_in_ux(i,3)~=0)
            g_x=out_index_x(i);
            g_y=out_index_y(i);
            if(point_material(g_x+1,g_y-1)==-1)
                R_D_x=d_d(g_x+1,g_y-1,1);
                R_D_y=d_d(g_x+1,g_y-1,2);
            else
                index_R_D=out_fic_index(g_x,g_y-1);
                for k=1:2
                    if(out_dir_index(index_R_D,k)==1)
                        R_D_x=out_fic_for_in_ux(index_R_D,k);
                        R_D_y=out_fic_for_in_uy(index_R_D,k);
                        break;
                    end
                end
            end
            if(point_material(g_x-1,g_y+1)==-1)
                L_U_x=d_d(g_x-1,g_y+1,1);
                L_U_y=d_d(g_x-1,g_y+1,2);
            else
                index_L_U=out_fic_index(g_x-1,g_y);
                for k=1:2
                    if(out_dir_index(index_L_U,k)==2)
                        L_U_x=out_fic_for_in_ux(index_L_U,k);
                        L_U_y=out_fic_for_in_uy(index_L_U,k);
                        break;
                    end
                end
            end
            out_fic_for_in_ux(i,3)=(4*out_fic_for_in_ux(i,1)+4*out_fic_for_in_ux(i,2)-4*d_d(g_x,g_y,1)+...
                (-L_U_x-R_D_x+d_d(g_x-1,g_y-1,1)))/3;
            out_fic_for_in_uy(i,3)=(4*out_fic_for_in_uy(i,1)+4*out_fic_for_in_uy(i,2)-4*d_d(g_x,g_y,2)+...
                (-L_U_y-R_D_y+d_d(g_x-1,g_y-1,2)))/3;
            %}
        end
    end
    %}
    %% new boundary nonlocal derivative computation
    for i=1:1:max_x
        for j=1:1:max_y
            t=point_material(i,j);
            if(t==-1)
                boun_point_type=out_fic_for_in(i,j);
                results=integer_displacement_gradient(i,j,t,boun_point_type,max_x,max_y);
            elseif(t==-2)
                results=integer_displacement_gradient(i,j,t,0,max_x,max_y);
            elseif(t==1)
                boun_point_type=in_fic_for_out(i,j);
                results=integer_displacement_gradient(i,j,t,boun_point_type,max_x,max_y);
            elseif(t==2)
                continue;
            end
            disp_gradient_x(i,j,1)=results(1);
            disp_gradient_x(i,j,2)=results(2);
            disp_gradient_y(i,j,1)=results(3);
            disp_gradient_y(i,j,2)=results(4);
        end
    end
    one_order_precision_disp_gradient(high_x,high_y,d_x,d_y);
    
    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        in_boun_1=in_cross_boun_mark(global_x,global_y,1);
        in_boun_2=in_cross_boun_mark(global_x,global_y,2);
        in_boun_disp_grad(i,in_boun_1,in_boun_2,global_x,global_y,size_out_fic_for_in);
    end
    
    frac_deri_x_disp_x_l=sym(zeros(high_x,high_y));
    frac_deri_x_disp_y_l=sym(zeros(high_x,high_y));
    frac_deri_x_disp_x_r=sym(zeros(high_x,high_y));
    frac_deri_x_disp_y_r=sym(zeros(high_x,high_y));
    frac_deri_y_disp_x_d=sym(zeros(high_x,high_y));
    frac_deri_y_disp_x_u=sym(zeros(high_x,high_y));
    frac_deri_y_disp_y_d=sym(zeros(high_x,high_y));
    frac_deri_y_disp_y_u=sym(zeros(high_x,high_y));
    global frac_deri_x_disp_x
    global frac_deri_x_disp_y
    global frac_deri_y_disp_x
    global frac_deri_y_disp_y
    frac_deri_x_disp_x=sym(zeros(high_x,high_y));
    frac_deri_x_disp_y=sym(zeros(high_x,high_y));
    frac_deri_y_disp_x=sym(zeros(high_x,high_y));
    frac_deri_y_disp_y=sym(zeros(high_x,high_y));
    % compute nonlocal derivative 
    for i=1:1:high_x
        for j=1:1:high_y
            if(point_material(i,j)<0)
                x_posi=(i-1)*d_x;
                y_posi=(j-1)*d_y;
                results=Caputo_derivative(left_info(i,j,:),i,j,1);
                lleft_x=results(1);
                lleft_y=results(2);
                results=Caputo_derivative(right_info(i,j,:),i,j,2);
                lright_x=results(1);
                lright_y=results(2);
                results=Caputo_derivative(down_info(i,j,:),i,j,3);
                ldown_x=results(1);
                ldown_y=results(2);
                results=Caputo_derivative(upper_info(i,j,:),i,j,4);
                lupper_x=results(1);
                lupper_y=results(2);
                if(left_info(i,j,2)==0)
                    frac_deri_x_disp_x_l(i,j)=lleft_x/GM_alpha;
                    frac_deri_x_disp_y_l(i,j)=lleft_y/GM_alpha;
                else
                    frac_deri_x_disp_x_l(i,j)=lleft_x*(x_posi-left_info(i,j,1))^(alpha-1);
                    frac_deri_x_disp_y_l(i,j)=lleft_y*(x_posi-left_info(i,j,1))^(alpha-1);
                end
                if(right_info(i,j,2)==0)
                    frac_deri_x_disp_x_r(i,j)=lright_x/GM_alpha;
                    frac_deri_x_disp_y_r(i,j)=lright_y/GM_alpha;
                else
                    frac_deri_x_disp_x_r(i,j)=lright_x*(right_info(i,j,1)-x_posi)^(alpha-1);
                    frac_deri_x_disp_y_r(i,j)=lright_y*(right_info(i,j,1)-x_posi)^(alpha-1);
                end
                if(down_info(i,j,2)==0)
                    frac_deri_y_disp_x_d(i,j)=ldown_x/GM_alpha;
                    frac_deri_y_disp_y_d(i,j)=ldown_y/GM_alpha;
                else
                    frac_deri_y_disp_x_d(i,j)=ldown_x*(y_posi-down_info(i,j,1))^(alpha-1);
                    frac_deri_y_disp_y_d(i,j)=ldown_y*(y_posi-down_info(i,j,1))^(alpha-1);
                end
                if(upper_info(i,j,2)==0)
                    frac_deri_y_disp_x_u(i,j)=lupper_x/GM_alpha;
                    frac_deri_y_disp_y_u(i,j)=lupper_y/GM_alpha;
                else
                    frac_deri_y_disp_x_u(i,j)=lupper_x*(upper_info(i,j,1)-y_posi)^(alpha-1);
                    frac_deri_y_disp_y_u(i,j)=lupper_y*(upper_info(i,j,1)-y_posi)^(alpha-1);
                end
                frac_deri_x_disp_x(i,j)=vpa((GM_alpha*(frac_deri_x_disp_x_l(i,j)+frac_deri_x_disp_x_r(i,j)))/2);
                frac_deri_x_disp_y(i,j)=vpa((GM_alpha*(frac_deri_x_disp_y_l(i,j)+frac_deri_x_disp_y_r(i,j)))/2);
                frac_deri_y_disp_x(i,j)=vpa((GM_alpha*(frac_deri_y_disp_x_d(i,j)+frac_deri_y_disp_x_u(i,j)))/2);
                frac_deri_y_disp_y(i,j)=vpa((GM_alpha*(frac_deri_y_disp_y_d(i,j)+frac_deri_y_disp_y_u(i,j)))/2);
            end
        end
    end    

    global sec_order_frac_deri_x_disp_x
    global sec_order_frac_deri_x_disp_y
    global sec_order_frac_deri_y_disp_x
    global sec_order_frac_deri_y_disp_y
    sec_order_frac_deri_x_disp_x=sym(zeros(high_x,high_y));
    sec_order_frac_deri_x_disp_y=sym(zeros(high_x,high_y));
    sec_order_frac_deri_y_disp_x=sym(zeros(high_x,high_y));
    sec_order_frac_deri_y_disp_y=sym(zeros(high_x,high_y));
    % compute nonlocal derivative 
    for i=1:1:high_x
        for j=1:1:high_y
            if(point_material(i,j)<0)
                x_posi=(i-1)*d_x;
                y_posi=(j-1)*d_y;
                results=sec_order_Caputo_derivative(left_info(i,j,:),i,j,1);
                lleft_x=results(1);
                lleft_y=results(2);
                results=sec_order_Caputo_derivative(right_info(i,j,:),i,j,2);
                lright_x=results(1);
                lright_y=results(2);
                results=sec_order_Caputo_derivative(down_info(i,j,:),i,j,3);
                ldown_x=results(1);
                ldown_y=results(2);
                results=sec_order_Caputo_derivative(upper_info(i,j,:),i,j,4);
                lupper_x=results(1);
                lupper_y=results(2);
                if(left_info(i,j,2)==0)
                    frac_deri_x_disp_x_l(i,j)=lleft_x/GM_alpha;
                    frac_deri_x_disp_y_l(i,j)=lleft_y/GM_alpha;
                else
                    frac_deri_x_disp_x_l(i,j)=lleft_x*(x_posi-left_info(i,j,1))^(alpha-1);
                    frac_deri_x_disp_y_l(i,j)=lleft_y*(x_posi-left_info(i,j,1))^(alpha-1);
                end
                if(right_info(i,j,2)==0)
                    frac_deri_x_disp_x_r(i,j)=lright_x/GM_alpha;
                    frac_deri_x_disp_y_r(i,j)=lright_y/GM_alpha;
                else
                    frac_deri_x_disp_x_r(i,j)=lright_x*(right_info(i,j,1)-x_posi)^(alpha-1);
                    frac_deri_x_disp_y_r(i,j)=lright_y*(right_info(i,j,1)-x_posi)^(alpha-1);
                end
                if(down_info(i,j,2)==0)
                    frac_deri_y_disp_x_d(i,j)=ldown_x/GM_alpha;
                    frac_deri_y_disp_y_d(i,j)=ldown_y/GM_alpha;
                else
                    frac_deri_y_disp_x_d(i,j)=ldown_x*(y_posi-down_info(i,j,1))^(alpha-1);
                    frac_deri_y_disp_y_d(i,j)=ldown_y*(y_posi-down_info(i,j,1))^(alpha-1);
                end
                if(upper_info(i,j,2)==0)
                    frac_deri_y_disp_x_u(i,j)=lupper_x/GM_alpha;
                    frac_deri_y_disp_y_u(i,j)=lupper_y/GM_alpha;
                else
                    frac_deri_y_disp_x_u(i,j)=lupper_x*(upper_info(i,j,1)-y_posi)^(alpha-1);
                    frac_deri_y_disp_y_u(i,j)=lupper_y*(upper_info(i,j,1)-y_posi)^(alpha-1);
                end
                sec_order_frac_deri_x_disp_x(i,j)=vpa((GM_alpha*(frac_deri_x_disp_x_l(i,j)+frac_deri_x_disp_x_r(i,j)))/2);
                sec_order_frac_deri_x_disp_y(i,j)=vpa((GM_alpha*(frac_deri_x_disp_y_l(i,j)+frac_deri_x_disp_y_r(i,j)))/2);
                sec_order_frac_deri_y_disp_x(i,j)=vpa((GM_alpha*(frac_deri_y_disp_x_d(i,j)+frac_deri_y_disp_x_u(i,j)))/2);
                sec_order_frac_deri_y_disp_y(i,j)=vpa((GM_alpha*(frac_deri_y_disp_y_d(i,j)+frac_deri_y_disp_y_u(i,j)))/2);
            end
        end
    end
    % recompute nonlocal derivative
    for i=1:1:size_out_fic_for_in
        %have problems here
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        for k=1:1:2
            s=out_dir_index(i,k);
            if(s~=0)
                result_l=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_left_info(i,k,:),i,k,global_x,global_y,1);
                result_r=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_right_info(i,k,:),i,k,global_x,global_y,2);
                result_d=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_down_info(i,k,:),i,k,global_x,global_y,3);
                result_u=boun_Caputo_derivative(boun_point_posi(i,k,:),boun_upper_info(i,k,:),i,k,global_x,global_y,4);
                if(boun_left_info(i,k,2)==0)
                    boun_frac_deri_x_disp_x_l(i,k)=result_l(1)/GM_alpha;
                    boun_frac_deri_x_disp_y_l(i,k)=result_l(2)/GM_alpha;
                else
                    boun_frac_deri_x_disp_x_l(i,k)=result_l(1)*(boun_point_posi(i,k,1)-boun_left_info(i,k,1))^(alpha-1);
                    boun_frac_deri_x_disp_y_l(i,k)=result_l(2)*(boun_point_posi(i,k,1)-boun_left_info(i,k,1))^(alpha-1);
                end
                if(boun_right_info(i,k,2)==0)
                    boun_frac_deri_x_disp_x_r(i,k)=result_r(1)/GM_alpha;
                    boun_frac_deri_x_disp_y_r(i,k)=result_r(2)/GM_alpha;
                else
                    boun_frac_deri_x_disp_x_r(i,k)=result_r(1)*(boun_right_info(i,k,1)-boun_point_posi(i,k,1))^(alpha-1);
                    boun_frac_deri_x_disp_y_r(i,k)=result_r(2)*(boun_right_info(i,k,1)-boun_point_posi(i,k,1))^(alpha-1);
                end
                if(boun_down_info(i,k,2)==0)
                    boun_frac_deri_y_disp_x_d(i,k)=result_d(1)/GM_alpha;
                    boun_frac_deri_y_disp_y_d(i,k)=result_d(2)/GM_alpha;
                else
                    boun_frac_deri_y_disp_x_d(i,k)=result_d(1)*(boun_point_posi(i,k,2)-boun_down_info(i,k,1))^(alpha-1);
                    boun_frac_deri_y_disp_y_d(i,k)=result_d(2)*(boun_point_posi(i,k,2)-boun_down_info(i,k,1))^(alpha-1);
                end
                if(boun_upper_info(i,k,2)==0)
                    boun_frac_deri_y_disp_x_u(i,k)=result_u(1)/GM_alpha;
                    boun_frac_deri_y_disp_y_u(i,k)=result_u(2)/GM_alpha;
                else
                    boun_frac_deri_y_disp_x_u(i,k)=result_u(1)*(boun_upper_info(i,k,1)-boun_point_posi(i,k,2))^(alpha-1);
                    boun_frac_deri_y_disp_y_u(i,k)=result_u(2)*(boun_upper_info(i,k,1)-boun_point_posi(i,k,2))^(alpha-1);
                end
            else
                continue;
            end
            boun_frac_deri_x_disp_x(i,k)=vpa((GM_alpha*(boun_frac_deri_x_disp_x_l(i,k)+boun_frac_deri_x_disp_x_r(i,k)))/2);
            boun_frac_deri_x_disp_y(i,k)=vpa((GM_alpha*(boun_frac_deri_x_disp_y_l(i,k)+boun_frac_deri_x_disp_y_r(i,k)))/2);
            boun_frac_deri_y_disp_x(i,k)=vpa((GM_alpha*(boun_frac_deri_y_disp_x_d(i,k)+boun_frac_deri_y_disp_x_u(i,k)))/2);
            boun_frac_deri_y_disp_y(i,k)=vpa((GM_alpha*(boun_frac_deri_y_disp_y_d(i,k)+boun_frac_deri_y_disp_y_u(i,k)))/2);
        end
    end
    
    global sym_boun_frac_deri_x_disp_x
    global sym_boun_frac_deri_x_disp_y
    global sym_boun_frac_deri_y_disp_x
    global sym_boun_frac_deri_y_disp_y
    sym_boun_frac_deri_x_disp_x=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_x_disp_y=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_x=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_y=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_x_disp_x_l=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_x_disp_y_l=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_x_d=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_y_d=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_x_disp_x_r=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_x_disp_y_r=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_x_u=sym(zeros(size_out_fic_for_in,2));
    sym_boun_frac_deri_y_disp_y_u=sym(zeros(size_out_fic_for_in,2));
    for i=1:1:size_out_fic_for_in
        %have problems here
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        for k=1:1:2
            s=out_dir_index(i,k);
            if(s~=0)
                result_l=new_sym_boun_Caputo_derivative(sym_boun_point_posi(i,k,:),sym_left_info(i,k,:),global_x,global_y,1);
                result_r=new_sym_boun_Caputo_derivative(sym_boun_point_posi(i,k,:),sym_right_info(i,k,:),global_x,global_y,2);
                result_d=new_sym_boun_Caputo_derivative(sym_boun_point_posi(i,k,:),sym_down_info(i,k,:),global_x,global_y,3);
                result_u=new_sym_boun_Caputo_derivative(sym_boun_point_posi(i,k,:),sym_upper_info(i,k,:),global_x,global_y,4);
                if(sym_left_info(i,k,2)==0)
                    sym_boun_frac_deri_x_disp_x_l(i,k)=result_l(1)/GM_alpha;
                    sym_boun_frac_deri_x_disp_y_l(i,k)=result_l(2)/GM_alpha;
                else
                    sym_boun_frac_deri_x_disp_x_l(i,k)=result_l(1)*(sym_boun_point_posi(i,k,1)-sym_left_info(i,k,1))^(alpha-1);
                    sym_boun_frac_deri_x_disp_y_l(i,k)=result_l(2)*(sym_boun_point_posi(i,k,1)-sym_left_info(i,k,1))^(alpha-1);
                end
                if(sym_right_info(i,k,2)==0)
                    sym_boun_frac_deri_x_disp_x_r(i,k)=result_r(1)/GM_alpha;
                    sym_boun_frac_deri_x_disp_y_r(i,k)=result_r(2)/GM_alpha;
                else
                    sym_boun_frac_deri_x_disp_x_r(i,k)=result_r(1)*(sym_right_info(i,k,1)-sym_boun_point_posi(i,k,1))^(alpha-1);
                    sym_boun_frac_deri_x_disp_y_r(i,k)=result_r(2)*(sym_right_info(i,k,1)-sym_boun_point_posi(i,k,1))^(alpha-1);
                end
                if(sym_down_info(i,k,2)==0)
                    sym_boun_frac_deri_y_disp_x_d(i,k)=result_d(1)/GM_alpha;
                    sym_boun_frac_deri_y_disp_y_d(i,k)=result_d(2)/GM_alpha;
                else
                    sym_boun_frac_deri_y_disp_x_d(i,k)=result_d(1)*(sym_boun_point_posi(i,k,2)-sym_down_info(i,k,1))^(alpha-1);
                    sym_boun_frac_deri_y_disp_y_d(i,k)=result_d(2)*(sym_boun_point_posi(i,k,2)-sym_down_info(i,k,1))^(alpha-1);
                end
                if(sym_upper_info(i,k,2)==0)
                    sym_boun_frac_deri_y_disp_x_u(i,k)=result_u(1)/GM_alpha;
                    sym_boun_frac_deri_y_disp_y_u(i,k)=result_u(2)/GM_alpha;
                else
                    sym_boun_frac_deri_y_disp_x_u(i,k)=result_u(1)*(sym_upper_info(i,k,1)-sym_boun_point_posi(i,k,2))^(alpha-1);
                    sym_boun_frac_deri_y_disp_y_u(i,k)=result_u(2)*(sym_upper_info(i,k,1)-sym_boun_point_posi(i,k,2))^(alpha-1);
                end
            else
                continue;
            end
            sym_boun_frac_deri_x_disp_x(i,k)=vpa((GM_alpha*(sym_boun_frac_deri_x_disp_x_l(i,k)+sym_boun_frac_deri_x_disp_x_r(i,k)))/2);
            sym_boun_frac_deri_x_disp_y(i,k)=vpa((GM_alpha*(sym_boun_frac_deri_x_disp_y_l(i,k)+sym_boun_frac_deri_x_disp_y_r(i,k)))/2);
            sym_boun_frac_deri_y_disp_x(i,k)=vpa((GM_alpha*(sym_boun_frac_deri_y_disp_x_d(i,k)+sym_boun_frac_deri_y_disp_x_u(i,k)))/2);
            sym_boun_frac_deri_y_disp_y(i,k)=vpa((GM_alpha*(sym_boun_frac_deri_y_disp_y_d(i,k)+sym_boun_frac_deri_y_disp_y_u(i,k)))/2);
        end
    end
    %}
    
    
    %% equilibrium fonfiguration
    global equilibrium
    equilibrium=sym(zeros(max_x,max_y,2));
    normal_point_equilibrium(max_x,max_y,miu_aluminum,lambda_aluminum,v_aluminum,miu_steel,lambda_steel,v_steel);
    out_boun_equilibrium(size_in_fic_for_out,miu_steel,lambda_steel,v_steel,d_x,d_y);
    in_boun_equilibrium(size_out_fic_for_in,miu_aluminum,lambda_aluminum,v_aluminum,d_x,d_y); 
    overall_boun_equilibrium(high_x,high_y,max_x,max_y,miu_aluminum,miu_steel,lambda_steel,d_x,d_y,size_in_fic_for_out,size_out_fic_for_in);

    %{
    % special treatment to type two boundary points
    index_1=1;
    global type_two_points
    type_two_points=zeros(1,3);
    for i=1:1:size_out_fic_for_in
        global_x=out_index_x(i);
        global_y=out_index_y(i);
        if(in_cross_boun_mark(global_x,global_y,1)==2)
            type_two_points(index_1,1)=global_x;
            type_two_points(index_1,2)=global_y;
            type_two_points(index_1,3)=i;
            index_1=index_1+1;
        end
    end
    size_type_two=size(type_two_points);
    for k=1:1:size_type_two
        i=type_two_points(k,1);
        j=type_two_points(k,2);
        k_t=type_two_points(k,3);
        disp_gradient_y(i+1,j,1)=(out_fic_for_in_ux(k_t,1)-d_d(i+1,j-1,1))/(2*d_y);
        disp_gradient_y(i+1,j,2)=(out_fic_for_in_uy(k_t,1)-d_d(i+1,j-1,2))/(2*d_y);
        disp_gradient_x(i,j+1,1)=(out_fic_for_in_ux(k_t,1)-d_d(i-1,j+1,1))/(2*d_x);
        disp_gradient_x(i,j+1,2)=(out_fic_for_in_uy(k_t,1)-d_d(i-1,j+1,2))/(2*d_x);
    end
    type_two_points_configuration(size_type_two,d_x,d_y,r);
    
    for k=1:1:size_type_two
        i=type_two_points(k,1);
        j=type_two_points(k,2);
        y_posi=(j-1)*d_y;
        results=sec_order_Caputo_derivative(down_info(i+1,j,:),i+1,j,3);
        ldown_x=results(1);
        ldown_y=results(2);
        results=sec_order_Caputo_derivative(upper_info(i+1,j,:),i+1,j,4);
        lupper_x=results(1);
        lupper_y=results(2);
        if(down_info(i+1,j,2)==0)
            frac_deri_y_disp_x_d(i+1,j)=ldown_x/GM_alpha;
            frac_deri_y_disp_y_d(i+1,j)=ldown_y/GM_alpha;
        else
            frac_deri_y_disp_x_d(i+1,j)=ldown_x*(y_posi-down_info(i+1,j,1))^(alpha-1);
            frac_deri_y_disp_y_d(i+1,j)=ldown_y*(y_posi-down_info(i+1,j,1))^(alpha-1);
        end
        if(upper_info(i+1,j,2)==0)
            frac_deri_y_disp_x_u(i+1,j)=lupper_x/GM_alpha;
            frac_deri_y_disp_y_u(i+1,j)=lupper_y/GM_alpha;
        else
            frac_deri_y_disp_x_u(i+1,j)=lupper_x*(upper_info(i+1,j,1)-y_posi)^(alpha-1);
            frac_deri_y_disp_y_u(i+1,j)=lupper_y*(upper_info(i+1,j,1)-y_posi)^(alpha-1);
        end
        sec_order_frac_deri_y_disp_x(i+1,j)=(GM_alpha*(frac_deri_y_disp_x_d(i+1,j)+frac_deri_y_disp_x_u(i+1,j)))/2;
        sec_order_frac_deri_y_disp_y(i+1,j)=(GM_alpha*(frac_deri_y_disp_y_d(i+1,j)+frac_deri_y_disp_y_u(i+1,j)))/2;
        
        x_posi=(i-1)*d_x;
        results=sec_order_Caputo_derivative(left_info(i,j+1,:),i,j+1,1);
        lleft_x=results(1);
        lleft_y=results(2);
        results=sec_order_Caputo_derivative(right_info(i,j+1,:),i,j+1,2);
        lright_x=results(1);
        lright_y=results(2);
        if(left_info(i,j+1,2)==0)
            frac_deri_x_disp_x_l(i,j+1)=lleft_x/GM_alpha;
            frac_deri_x_disp_y_l(i,j+1)=lleft_y/GM_alpha;
        else
            frac_deri_x_disp_x_l(i,j+1)=lleft_x*(x_posi-left_info(i,j+1,1))^(alpha-1);
            frac_deri_x_disp_y_l(i,j+1)=lleft_y*(x_posi-left_info(i,j+1,1))^(alpha-1);
        end
        if(right_info(i,j+1,2)==0)
            frac_deri_x_disp_x_r(i,j+1)=lright_x/GM_alpha;
            frac_deri_x_disp_y_r(i,j+1)=lright_y/GM_alpha;
        else
            frac_deri_x_disp_x_r(i,j+1)=lright_x*(right_info(i,j+1,1)-x_posi)^(alpha-1);
            frac_deri_x_disp_y_r(i,j+1)=lright_y*(right_info(i,j+1,1)-x_posi)^(alpha-1);
        end
        sec_order_frac_deri_x_disp_x(i,j+1)=(GM_alpha*(frac_deri_x_disp_x_l(i,j+1)+frac_deri_x_disp_x_r(i,j+1)))/2;
        sec_order_frac_deri_x_disp_y(i,j+1)=(GM_alpha*(frac_deri_x_disp_y_l(i,j+1)+frac_deri_x_disp_y_r(i,j+1)))/2;
    end
    
    for k=1:1:size_type_two
        m=type_two_points(k,1);
        n=type_two_points(k,2);
        equilibrium(m,n,1)=vpa((miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
            (frac_deri_x_disp_x(m,n)-frac_deri_x_disp_x(m-1,n))/(d_x)+...
            (1-2*v_aluminum)*(frac_deri_y_disp_x(m,n)-frac_deri_y_disp_x(m,n-1))/(d_y)+...
            (1-2*v_aluminum)*(sec_order_frac_deri_x_disp_y(m,n+1)-sec_order_frac_deri_x_disp_y(m,n-1))/(2*d_y)+...
            2*v_aluminum*(sec_order_frac_deri_y_disp_y(m+1,n)-sec_order_frac_deri_y_disp_y(m-1,n))/(2*d_x)));
        equilibrium(m,n,2)=vpa((miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
            (frac_deri_y_disp_y(m,n)-frac_deri_y_disp_y(m,n-1))/(d_y)+...
            (1-2*v_aluminum)*(frac_deri_x_disp_y(m,n)-frac_deri_x_disp_y(m-1,n))/(d_x)+...
            (1-2*v_aluminum)*(sec_order_frac_deri_y_disp_x(m+1,n)-sec_order_frac_deri_y_disp_x(m-1,n))/(2*d_x)+...
            2*v_aluminum*(sec_order_frac_deri_x_disp_x(m,n+1)-sec_order_frac_deri_x_disp_x(m,n-1))/(2*d_y)));
        
    end
    %}
    
    %% final compuation configuration
    %{
    symbolic_variable_cell=cell(max_x,max_y,2);
	variable_subscript_cell=cell(max_x,max_y,2);
	for i=1:1:max_x
		for j=1:1:max_y
			for k=1:1:2
				symbolic_variables=symvar(equilibrium(i,j,k));
				symbolic_variable_cell(i,j,k)={symbolic_variables};
				number_of_variable=length(symbolic_variables);
				subscript_number=zeros(number_of_variable,3);
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
						subscript_number(p,3)=1;
					else
						subscript_number(p,3)=2;
					end
					string_1=str_variable(place(1)+1:1:place(2)-1);
					string_2=str_variable(place(2)+1:1:length_of_string);
					subscript_number(p,1)=str2double(string_1);
					subscript_number(p,2)=str2double(string_2);  
				end
				variable_subscript_cell(i,j,k)={subscript_number};  
			end
		end
	end
	
	coefficient_displacement=zeros(max_x,max_y,2);
	%compute corresponding coefficient of each displacement Ux_i_j and Uy_i_j 

	for i=1:1:max_x
		for j=1:1:max_y
			if((i==1)&&(j==1))
				coefficient_displacement(i,j,1)=0;
				coefficient_displacement(i,j,2)=0;
				
			elseif((i==1)&&(j~=1))
				coefficient_displacement(i,j,1)=0;
				coeff=coeffs(equilibrium(i,j,2),Uy(i,j));
				if(size(coeff,2)==1)
					coefficient_displacement(i,j,2)=coeff(1);
				elseif(size(coeff,2)==2)
					coefficient_displacement(i,j,2)=coeff(2);
				end

			elseif((i~=1)&&(j==1))
				coeff=coeffs(equilibrium(i,j,1),Ux(i,j));
				if(size(coeff,2)==1)
					coefficient_displacement(i,j,1)=coeff(1);
				elseif(size(coeff,2)==2)
					coefficient_displacement(i,j,1)=coeff(2);
				end
				coefficient_displacement(i,j,2)=0;
				
			else
				coeff=coeffs(equilibrium(i,j,1),Ux(i,j));
				if(size(coeff,2)==1)
					coefficient_displacement(i,j,1)=coeff(1);
				elseif(size(coeff,2)==2)
					coefficient_displacement(i,j,1)=coeff(2);
				end
				
				coeff=coeffs(equilibrium(i,j,2),Uy(i,j));
				if(size(coeff,2)==1)
					coefficient_displacement(i,j,2)=coeff(1);
				elseif(size(coeff,2)==2)
					coefficient_displacement(i,j,2)=coeff(2);
				end
			end
		end
	end
    abc=zeros(high_x,high_y,2);
    for i=1:1:high_x
        for j=1:1:high_y
            abc(i,j,1)=vpa(coefficient_displacement(i,j,1)-coefficient_displacement(j,i,2));
            abc(i,j,2)=vpa(coefficient_displacement(i,j,2)-coefficient_displacement(j,i,1));
        end
    end
    
	coefficient_b=zeros(max_x,max_y,2);
    for i=1:1:max_x
        for j=1:1:max_y
            if(i==max_x)
                coefficient_b(i,j,2)=0;%1000N force added on the boundary
                coefficient_b(i,j,1)=1;
                coefficient_b1(i+max_y*(j-1))=1;
            else
                coefficient_b(i,j,1)=0;
                coefficient_b(i,j,2)=0;
            end
        end
    end
	%iterations
	accuracy=5e-16;
    accuracy=1e-7;
    global displacement
	displacement=zeros(max_x,max_y,2);
	updated_displacement=zeros(max_x,max_y,2);
	%method to parse variables in the equilibrium before enter into the
	%iteration process
    
	variable_coefficient_cell=cell(max_x,max_y,2);
    for i=1:1:max_x
        for j=1:1:max_y
            for k=1:1:2
                variables=symbolic_variable_cell{i,j,k};
                if(isempty(variables))
                    continue;
                else
                    coeff_variables=zeros(length(variables),1);
                    %get coefficient of each variables in each equilibrium
                    variables=symvar(equilibrium(i,j,k));
                    coeff=flip(coeffs(equilibrium(i,j,k),variables));
                    for m=1:1:length(variables)
                        coeff_variables(m)=coeff(m);
                    end
                    variable_coefficient_cell{i,j,k}=coeff_variables;
                end
            end
        end
    end
    
    tic;
    %final computation
    iteration_accuracy=10;
    while(iteration_accuracy>accuracy)
        T=0;
        for i=1:1:max_x
            for j=1:1:max_y
                for k=1:1:2
                    if(coefficient_displacement(i,j,k)~=0)
                        %select out those variables that participate in the
                        %computation
                        variables=symbolic_variable_cell{i,j,k};
                        subscript_variables=variable_subscript_cell{i,j,k};
                        coeff_variables=variable_coefficient_cell{i,j,k};
                        residue=0;
                        for m=1:1:size(variables,2)
                            residue=residue+displacement(subscript_variables(m,1),subscript_variables(m,2),subscript_variables(m,3))*coeff_variables(m);
                        end       
                        updated_displacement(i,j,k)=displacement(i,j,k)-1.1*(residue-coefficient_b(i,j,k))/coefficient_displacement(i,j,k);
                        %if(displacement(i,j,k)~=0)
                        T=max(abs((updated_displacement(i,j,k)-displacement(i,j,k))),T);
                        %end
                        displacement(i,j,k)=updated_displacement(i,j,k);
                    end
                end
            end
        end
        iteration_accuracy=T;
    end
    toc;
%}
    %% New iteration method: Weighted Jacobian iteration
    var_coeff_matrix=sparse(max_x*max_y*2, max_x*max_y*2);
    for i=1:1:max_x
        for j=1:1:max_y
            for k=1:1:2
                [var_coeffs, vars]=coeffs(equilibrium(i,j,k));
                for p=1:1:length(vars)
                    str_variable=char(vars(p));
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
                    string_1=str_variable(place(1)+1:1:place(2)-1);
					string_2=str_variable(place(2)+1:1:length_of_string);
					sub_num_1=str2double(string_1);
					sub_num_2=str2double(string_2);  
                    if(str_variable(2)=='x')% 1 stands for Ux; 2 stands for Uy;
						var_coeff_matrix(i+max_y*(j-1)+max_x*max_y*(k-1),sub_num_1+max_y*(sub_num_2-1))=var_coeffs(p);
					else
						var_coeff_matrix(i+max_y*(j-1)+max_x*max_y*(k-1),sub_num_1+max_y*(sub_num_2-1)+max_x*max_y)=var_coeffs(p);
                    end
                end
            end
        end
    end
    coefficient_b1=zeros(max_x*max_y*2,1);
    for i=1:1:max_x
        for j=1:1:max_y
            if(i==max_x)
                coefficient_b1(i+max_y*(j-1))=1;
            end
        end
    end
    
    Rela_residue=1e-4;
    
    %w=0.9;
    D=diag(var_coeff_matrix);
    L=-tril(var_coeff_matrix)+diag(D);
    U=diag(D)-L-var_coeff_matrix;
    %D=D/w;
    inv_D=zeros(size(D));
    for i=1:length(D)
        if(D(i)~=0)
            inv_D(i)=1/D(i);
        end
    end
    inv_D=diag(inv_D);
    D=diag(D);
    U_old=zeros(max_x*max_y*2,1);
    U_new=zeros(max_x*max_y*2,1);
    V=D-L;
    norm_b1=norm(coefficient_b1);
    tic;
    while(1)
        %U_new=inv_D*((D-var_coeff_matrix)*U_old+coefficient_b1);
        b=U*U_old+coefficient_b1;
        if(V(1,1)==0)
            U_new(1)=0;
        else
            U_new(1,1)=b(1)/V(1,1);
        end
        for i=2:length(coefficient_b1)
            if(V(i,i)==0)
                U_new(i)=0;
            else
                U_new(i)=(b(i)-V(i,1:i-1)*U_new(1:i-1))/V(i,i);
            end
        end
        residue=norm(coefficient_b1-var_coeff_matrix*U_new)/norm_b1;
        U_old=U_new;
        if(residue<Rela_residue)
            break;
        end
        if(norm(U_new)>1e8)
            break;
        end
    end
    toc;
    %% post process
    global stress_x
    global stress_y
    global stress_xy
    global x_position
    global y_position
    extra_y=2*ceil(size_out_fic_for_in/max_x)+2*ceil(size_in_fic_for_out/max_x);
    x_position=zeros(max_x,max_y+extra_y);
    y_position=zeros(max_x,max_y+extra_y);
    stress_x=zeros(max_x,max_y+extra_y);
    stress_y=zeros(max_x,max_y+extra_y);
    stress_xy=zeros(max_x,max_y+extra_y);
    for i=1:1:size_in_fic_for_out
        global_x=in_index_x(i);
        global_y=in_index_y(i);
        out_boun_1=out_cross_boun_mark(global_x,global_y,1);
        out_boun_2=out_cross_boun_mark(global_x,global_y,2);
        out_boun_disp_grad(i,out_boun_1,out_boun_2,global_x,global_y,size_in_fic_for_out);
    end
    out_deri_x_disp_x=vpa(out_deri_x_disp_x);
    out_deri_x_disp_y=vpa(out_deri_x_disp_y);
    out_deri_y_disp_x=vpa(out_deri_y_disp_x);
    out_deri_y_disp_y=vpa(out_deri_y_disp_y);
    stress_compute(r,d_x,d_y,high_x,high_y,max_x,max_y,miu_aluminum,lambda_aluminum,miu_steel,lambda_steel,size_out_fic_for_in,size_in_fic_for_out);
    x_posi_col=x_position(:);
    y_posi_col=y_position(:);
    x_posi_col=(x_posi_col-1)*d_x;
    y_posi_col=(y_posi_col-1)*d_y;
    stress_x_col=stress_x(:);
    stress_y_col=stress_y(:);
    stress_xy_col=stress_xy(:);
    total_disp=sqrt(displacement(:,:,1).^2+displacement(:,:,2).^2);
    x_posi=zeros(max_x,max_y);
    y_posi=zeros(max_x,max_y);
    for i=1:1:max_x
        x_posi(i,:)=(i-1)*d_x;
        y_posi(:,i)=(i-1)*d_y;
    end
    save 'new-0.8-4-0.05.mat';