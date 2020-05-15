function[]=quarter_circle_inclusion()
    %material properties
    v_steel=0.28883542806974205;
    miu_steel=8.237068454836098e10;
    lambda_steel=1.1266838804655762e11;
    %v_aluminum=0.28883542806974205;
    %miu_aluminum=8.237068454836098e10;
    %lambda_aluminum=1.1266838804655762e11;
    v_aluminum=0.33124956931271565;
    lambda_aluminum=5.097269526934259e10;
    miu_aluminum=2.596732215483388e10;
    %ggeometry configuration
    l_x=1.5;
    l_y=1.5;
    r=0.925;
    d_x=0.05;
    d_y=0.05;
    center_x=0;
    center_y=0;
    max_x=l_x/d_x+1;
    max_y=l_y/d_y+1;
    
    d_d=sym(zeros(max_x,max_y,2));
    %discretization grid
    Ux=sym('Ux_',[l_x/d_x+1 l_y/d_y+1]);
    Uy=sym('Uy_',[l_x/d_x+1 l_y/d_y+1]);
    for i=1:1:l_x/d_x+1
        for j=1:1:l_y/d_y+1
            d_d(i,j,1)=Ux(i,j);
            d_d(i,j,2)=Uy(i,j);
        end
    end
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
    
	in_cross_boun_mark=zeros(high_x,high_y,2);
	out_cross_boun_mark=zeros(high_x,high_y,2);
	
	% for outer boundary points
	%   type(1):
	%   x .
	%   x .    (1)
	%   x .
	%
	%   . . .
	%   x x x  (2)
	% type(2):
	%        * .  (1)
	%        * *   
	% type(3):
	%        . .  (1)
	%        * .
	% type(4):
	%        (1)
	%        . . . 
	%        * * .    
	%             . .
	%             * . (2)
	%             * .
	% type(5):
	%        (1)
	%        x . . 
	%        x x x    
	%             x .
	%             x . (2)
	%             x x
	
	
	% for inner boundary points
	% type(1)
	%   . x
	%   . x    (1)
	%   . x
	%   x x x
	%   . . .  (2)
	% type(2):
	%        . .  (1)
	%        x .   
	% type(3):
	%        x x  (1)
	%        . x
	% type(4):
	%        (1)
	%        x x x 
	%        . . x    
	%             x x
	%             . x (2)
	%             . x
	% type(5):
	%        (1)
	%        . x x 
	%        . . .    
	%             . x
	%             . x (2)
	%             . .
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
                in_fic_for_out_ux(index,3)=sym('in_fic_ux_'+string(i)+'_'+string(j)+'_LB');
                in_fic_for_out_uy(index,3)=sym('in_fic_uy_'+string(i)+'_'+string(j)+'_LB');
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
                out_fic_for_in_ux(index_1,3)=sym('out_fic_ux_'+string(i)+'_'+string(j)+'_RT');
                out_fic_for_in_uy(index_1,3)=sym('out_fic_uy_'+string(i)+'_'+string(j)+'_RT');
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
	
	in_fic_index=zeros(high_x,high_y);
	out_fic_index=zeros(high_x,high_y);
	
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
	
	size_in_fic_for_out=length(in_fic_for_out_ux);
	size_out_fic_for_in=length(out_fic_for_in_ux);
   
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
	%normal direction of the surface
	in_sur_nor_dir_x=zeros(size_in_fic_for_out,3);
	in_sur_nor_dir_y=zeros(size_in_fic_for_out,3);
	out_sur_nor_dir_x=zeros(size_out_fic_for_in,3);
	out_sur_nor_dir_y=zeros(size_out_fic_for_in,3);
    
    for i=1:1:size_in_fic_for_out
        if(in_fic_for_out_ux(i,2)==0)
            str_variable=char(in_fic_for_out_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1)=='_')
                if(str_variable(length_str)=='L')
                    dis_y=(in_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    in_sur_nor_dir_x(i,1)=dis_x-center_x;
                    in_sur_nor_dir_y(i,1)=dis_y;
                    %ux situation
                    x_disp=[in_fic_for_out_ux(i,1);d_d(in_index_x(i),in_index_y(i),1);d_d(in_index_x(i)+1,in_index_y(i),1)];
                    syms a b c;
                    eqn(1)=a*((in_index_x(i)-2)*d_x)^2+b*((in_index_x(i)-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((in_index_x(i)-1)*d_x)^2+b*((in_index_x(i)-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((in_index_x(i))*d_x)^2+b*((in_index_x(i))*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
                    out_deri_x_disp_x(i,1)=2*a*dis_x+b;
                    %uy situation
                    y_disp=[in_fic_for_out_uy(i,1);d_d(in_index_x(i),in_index_y(i),2);d_d(in_index_x(i)+1,in_index_y(i),2)];
                    syms a b c;
                    eqn(1)=a*((in_index_x(i)-2)*d_x)^2+b*((in_index_x(i)-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((in_index_x(i)-1)*d_x)^2+b*((in_index_x(i)-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((in_index_x(i))*d_x)^2+b*((in_index_x(i))*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
                    out_deri_x_disp_y(i,1)=2*a*dis_x+b;  
                end
                if(str_variable(length_str)=='B')
                    dis_x=(in_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y;
                    in_sur_nor_dir_x(i,1)=dis_x;
                    in_sur_nor_dir_y(i,1)=dis_y-center_y;
                    x_disp=[in_fic_for_out_ux(i,1);d_d(in_index_x(i),in_index_y(i),1);d_d(in_index_x(i),in_index_y(i)+1,1)];
                    syms a b c;
                    eqn(1)=a*((in_index_y(i)-2)*d_y)^2+b*((in_index_y(i)-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((in_index_y(i)-1)*d_y)^2+b*((in_index_y(i)-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((in_index_y(i))*d_y)^2+b*((in_index_y(i))*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
                    out_deri_y_disp_x(i,1)=2*a*dis_y+b;
                    %uy situation
                    y_disp=[in_fic_for_out_uy(i,1);d_d(in_index_x(i),in_index_y(i),2);d_d(in_index_x(i),in_index_y(i)+1,2)];
                    syms a b c;
                    eqn(1)=a*((in_index_y(i)-2)*d_y)^2+b*((in_index_y(i)-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((in_index_y(i)-1)*d_y)^2+b*((in_index_y(i)-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((in_index_y(i))*d_y)^2+b*((in_index_y(i))*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
                    out_deri_y_disp_y(i,1)=2*a*dis_y+b;
                end
            else
                continue;
            end			
        else
            for k=1:1:2
                str_variable=char(in_fic_for_out_ux(i,k));
                length_str=length(str_variable);
                if(str_variable(length_str)=='L')
                    dis_y=(in_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    in_sur_nor_dir_x(i,k)=dis_x-center_x;
                    in_sur_nor_dir_y(i,k)=dis_y;
                    %ux situation
                    x_disp=[in_fic_for_out_ux(i,k);d_d(in_index_x(i),in_index_y(i),1);d_d(in_index_x(i)+1,in_index_y(i),1)];
                    syms a b c;
                    eqn(1)=a*((in_index_x(i)-2)*d_x)^2+b*((in_index_x(i)-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((in_index_x(i)-1)*d_x)^2+b*((in_index_x(i)-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((in_index_x(i))*d_x)^2+b*((in_index_x(i))*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_x(i,k)=a*dis_x^2+b*dis_x+c;
                    out_deri_x_disp_x(i,k)=2*a*dis_x+b;
                    %uy situation
                    y_disp=[in_fic_for_out_uy(i,k);d_d(in_index_x(i),in_index_y(i),2);d_d(in_index_x(i)+1,in_index_y(i),2)];
                    syms a b c;
                    eqn(1)=a*((in_index_x(i)-2)*d_x)^2+b*((in_index_x(i)-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((in_index_x(i)-1)*d_x)^2+b*((in_index_x(i)-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((in_index_x(i))*d_x)^2+b*((in_index_x(i))*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_y(i,k)=a*dis_x^2+b*dis_x+c;
                    out_deri_x_disp_y(i,k)=2*a*dis_x+b;
                end
                if(str_variable(length_str)=='B')
                    dis_x=(in_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y;
                    in_sur_nor_dir_x(i,k)=dis_x;
                    in_sur_nor_dir_y(i,k)=dis_y-center_y;
                    x_disp=[in_fic_for_out_ux(i,k);d_d(in_index_x(i),in_index_y(i),1);d_d(in_index_x(i),in_index_y(i)+1,1)];
                    syms a b c;
                    eqn(1)=a*((in_index_y(i)-2)*d_y)^2+b*((in_index_y(i)-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((in_index_y(i)-1)*d_y)^2+b*((in_index_y(i)-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((in_index_y(i))*d_y)^2+b*((in_index_y(i))*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_x(i,k)=a*dis_y^2+b*dis_y+c;
                    out_deri_y_disp_x(i,k)=2*a*dis_y+b;
                    %uy situation
                    y_disp=[in_fic_for_out_uy(i,k);d_d(in_index_x(i),in_index_y(i),2);d_d(in_index_x(i),in_index_y(i)+1,2)];
                    syms a b c;
                    eqn(1)=a*((in_index_y(i)-2)*d_y)^2+b*((in_index_y(i)-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((in_index_y(i)-1)*d_y)^2+b*((in_index_y(i)-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((in_index_y(i))*d_y)^2+b*((in_index_y(i))*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_disp_y(i,k)=a*dis_y^2+b*dis_y+c;
                    out_deri_y_disp_y(i,k)=2*a*dis_y+b;
                end
            end
        end
    end
	
    for i=1:1:size_in_fic_for_out
        if(in_fic_for_out_ux(i,2)~=0)
            index_x=in_index_x(i);
            index_y=in_index_y(i);
            in_fic_for_out_ux(i,3)=in_fic_for_out_ux(i,1)+in_fic_for_out_ux(i,2)-d_d(index_x,index_y,1);
            in_fic_for_out_uy(i,3)=in_fic_for_out_uy(i,1)+in_fic_for_out_uy(i,2)-d_d(index_x,index_y,2);
        end
    end
    
	for i=1:1:size_out_fic_for_in
		if(out_fic_for_in_ux(i,2)==0)
			str_variable=char(out_fic_for_in_ux(i,1));
			length_str=length(str_variable);
			if(str_variable(length_str-1)=='_')
				if(str_variable(length_str)=='R')
					dis_y=(out_index_y(i)-1)*d_y-center_y;
					dis_x=sqrt(r^2-dis_y^2)+center_x;
					out_sur_nor_dir_x(i,1)=dis_x-center_x;
					out_sur_nor_dir_y(i,1)=dis_y;
					%ux direction
					x_disp=[d_d(out_index_x(i)-1,out_index_y(i),1);d_d(out_index_x(i),out_index_y(i),1);out_fic_for_in_ux(i,1)];
					syms a b c;
					eqn(1)=a*((out_index_x(i)-2)*d_x)^2+b*((out_index_x(i)-2)*d_x)+c==x_disp(1);
					eqn(2)=a*((out_index_x(i)-1)*d_x)^2+b*((out_index_x(i)-1)*d_x)+c==x_disp(2);
					eqn(3)=a*((out_index_x(i))*d_x)^2+b*((out_index_x(i))*d_x)+c==x_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_x(i,1)=a*dis_x^2+b*dis_x+c;
					in_deri_x_disp_x(i,1)=2*a*dis_x+b;
					%uy direction
					y_disp=[d_d(out_index_x(i)-1,out_index_y(i),2);d_d(out_index_x(i),out_index_y(i),2);out_fic_for_in_uy(i,1)];
					syms a b c;
					eqn(1)=a*((out_index_x(i)-2)*d_x)^2+b*((out_index_x(i)-2)*d_x)+c==y_disp(1);
					eqn(2)=a*((out_index_x(i)-1)*d_x)^2+b*((out_index_x(i)-1)*d_x)+c==y_disp(2);
					eqn(3)=a*((out_index_x(i))*d_x)^2+b*((out_index_x(i))*d_x)+c==y_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_y(i,1)=a*dis_x^2+b*dis_x+c;
					in_deri_x_disp_y(i,1)=2*a*dis_x+b;
				end
				if(str_variable(length_str)=='T')
					dis_x=(out_index_x(i)-1)*d_x-center_x;
					dis_y=sqrt(r^2-dis_x^2)+center_y;
					out_sur_nor_dir_x(i,1)=dis_x;
					out_sur_nor_dir_y(i,1)=dis_y-center_y;
					%ux direction
					x_disp=[d_d(out_index_x(i),out_index_y(i)-1,1);d_d(out_index_x(i),out_index_y(i),1);out_fic_for_in_ux(i,1)];
					syms a b c;
					eqn(1)=a*((out_index_y(i)-2)*d_y)^2+b*((out_index_y(i)-2)*d_y)+c==x_disp(1);
					eqn(2)=a*((out_index_y(i)-1)*d_y)^2+b*((out_index_y(i)-1)*d_y)+c==x_disp(2);
					eqn(3)=a*((out_index_y(i))*d_y)^2+b*((out_index_y(i))*d_y)+c==x_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_x(i,1)=a*dis_y^2+b*dis_y+c;
					in_deri_y_disp_x(i,1)=2*a*dis_y+b;
					%uy direction
					y_disp=[d_d(out_index_x(i),out_index_y(i)-1,2);d_d(out_index_x(i),out_index_y(i),2);out_fic_for_in_uy(i,1)];
					syms a b c;
					eqn(1)=a*((out_index_y(i)-2)*d_y)^2+b*((out_index_y(i)-2)*d_y)+c==y_disp(1);
					eqn(2)=a*((out_index_y(i)-1)*d_y)^2+b*((out_index_y(i)-1)*d_y)+c==y_disp(2);
					eqn(3)=a*((out_index_y(i))*d_y)^2+b*((out_index_y(i))*d_y)+c==y_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_y(i,1)=a*dis_y^2+b*dis_y+c;
					in_deri_y_disp_y(i,1)=2*a*dis_y+b;
				end
			else
				continue;
			end
		else
			for k=1:1:2
				str_variable=char(out_fic_for_in_ux(i,k));
				length_str=length(str_variable);
				if(str_variable(length_str)=='R')
					dis_y=(out_index_y(i)-1)*d_y-center_y;
					dis_x=sqrt(r^2-dis_y^2)+center_x;
					out_sur_nor_dir_x(i,k)=dis_x-center_x;
					out_sur_nor_dir_y(i,k)=dis_y;
					%ux direction
					x_disp=[d_d(out_index_x(i)-1,out_index_y(i),1);d_d(out_index_x(i),out_index_y(i),1);out_fic_for_in_ux(i,k)];
					syms a b c;
					eqn(1)=a*((out_index_x(i)-2)*d_x)^2+b*((out_index_x(i)-2)*d_x)+c==x_disp(1);
					eqn(2)=a*((out_index_x(i)-1)*d_x)^2+b*((out_index_x(i)-1)*d_x)+c==x_disp(2);
					eqn(3)=a*((out_index_x(i))*d_x)^2+b*((out_index_x(i))*d_x)+c==x_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_x(i,k)=a*dis_x^2+b*dis_x+c;
					in_deri_x_disp_x(i,k)=2*a*dis_x+b;
					%uy direction
					y_disp=[d_d(out_index_x(i)-1,out_index_y(i),2);d_d(out_index_x(i),out_index_y(i),2);out_fic_for_in_uy(i,k)];
					syms a b c;
					eqn(1)=a*((out_index_x(i)-2)*d_x)^2+b*((out_index_x(i)-2)*d_x)+c==y_disp(1);
					eqn(2)=a*((out_index_x(i)-1)*d_x)^2+b*((out_index_x(i)-1)*d_x)+c==y_disp(2);
					eqn(3)=a*((out_index_x(i))*d_x)^2+b*((out_index_x(i))*d_x)+c==y_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_y(i,k)=a*dis_x^2+b*dis_x+c;
					in_deri_x_disp_y(i,k)=2*a*dis_x+b;
				end
				if(str_variable(length_str)=='T')
					dis_x=(out_index_x(i)-1)*d_x-center_x;
					dis_y=sqrt(r^2-dis_x^2)+center_y;
					out_sur_nor_dir_x(i,k)=dis_x;
					out_sur_nor_dir_y(i,k)=dis_y-center_y;
					%ux direction
					x_disp=[d_d(out_index_x(i),out_index_y(i)-1,1);d_d(out_index_x(i),out_index_y(i),1);out_fic_for_in_ux(i,k)];
					syms a b c;
					eqn(1)=a*((out_index_y(i)-2)*d_y)^2+b*((out_index_y(i)-2)*d_y)+c==x_disp(1);
					eqn(2)=a*((out_index_y(i)-1)*d_y)^2+b*((out_index_y(i)-1)*d_y)+c==x_disp(2);
					eqn(3)=a*((out_index_y(i))*d_y)^2+b*((out_index_y(i))*d_y)+c==x_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_x(i,k)=a*dis_y^2+b*dis_y+c;
					in_deri_y_disp_x(i,k)=2*a*dis_y+b;
					%uy direction
					y_disp=[d_d(out_index_x(i),out_index_y(i)-1,2);d_d(out_index_x(i),out_index_y(i),2);out_fic_for_in_uy(i,k)];
					syms a b c;
					eqn(1)=a*((out_index_y(i)-2)*d_y)^2+b*((out_index_y(i)-2)*d_y)+c==y_disp(1);
					eqn(2)=a*((out_index_y(i)-1)*d_y)^2+b*((out_index_y(i)-1)*d_y)+c==y_disp(2);
					eqn(3)=a*((out_index_y(i))*d_y)^2+b*((out_index_y(i))*d_y)+c==y_disp(3);
					[a,b,c]=solve(eqn,[a,b,c]);
					in_disp_y(i,k)=a*dis_y^2+b*dis_y+c;
					in_deri_y_disp_y(i,k)=2*a*dis_y+b;
				end
			end
		end
	end

	for i=1:1:length(out_fic_for_in_ux)
		if(out_fic_for_in_ux(i,2)~=0)
			index_x=out_index_x(i);
			index_y=out_index_y(i);
			out_fic_for_in_ux(i,3)=out_fic_for_in_ux(i,1)+out_fic_for_in_ux(i,2)-d_d(index_x,index_y,1);
			out_fic_for_in_uy(i,3)=out_fic_for_in_uy(i,1)+out_fic_for_in_uy(i,2)-d_d(index_x,index_y,2);
		end
	end

	out_sur_nor_dir_x=vpa(out_sur_nor_dir_x);
	out_sur_nor_dir_y=vpa(out_sur_nor_dir_y);

	%boundary points
	x_global=in_index_x(1);
	y_global=in_index_y(1);
	x_local=x_global;
	y_local=y_global;
	dis_x=(in_index_x(1)-1)*d_x-center_x;
	dis_y=sqrt(r^2-dis_x^2)+center_y;
	%right2 situation
	j=in_fic_index(x_local+2,y_local);
    for k=1:1:2
        if(in_dir_index(j,k)==2)
            in_fic_right2_x=in_fic_for_out_ux(j,k);
            in_fic_right2_y=in_fic_for_out_uy(j,k);
            break;
        end
    end
	%compute right2 out_disp_x
	x_disp=[in_fic_right2_x;d_d(x_global+2,y_global,1);d_d(x_global+2,y_global+1,1)];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_right2_disp_x=a*dis_y^2+b*dis_y+c;
	%compute right2 out_disp_y
	y_disp=[in_fic_right2_y;d_d(x_global+2,y_global,2);d_d(x_global+2,y_global+1,2)];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_right2_disp_y=a*dis_y^2+b*dis_y+c;
	%right situation
	j=in_fic_index(x_local+1,y_local);
    for k=1:1:2
        if(in_dir_index(j,k)==2)
            in_fic_right_x=in_fic_for_out_ux(j,k);
            in_fic_right_y=in_fic_for_out_uy(j,k);
            break;
        end
    end
    
	%compute right out_disp_x
	x_disp=[in_fic_right_x;d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_right_disp_x=a*dis_y^2+b*dis_y+c;
	%compute right out_disp_y
	y_disp=[in_fic_right_y;d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_right_disp_y=a*dis_y^2+b*dis_y+c;
	out_deri_x_disp_x(1,1)=(-3*out_disp_x(1,1)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
	out_deri_x_disp_y(1,1)=(-3*out_disp_y(1,1)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
	
	
	%boundary points
	x_global=in_index_x(size_in_fic_for_out);
	y_global=in_index_y(size_in_fic_for_out);
	x_local=x_global;
	y_local=y_global;
	dis_y=(in_index_y(size_in_fic_for_out)-1)*d_y-center_y;
	dis_x=sqrt(r^2-dis_y^2)+center_x;
	%the index of in_fic on outer point (x,y+1)
	%the corresponding in_fic_point is on the left of (x,y+1)
	%upper situation
	j=in_fic_index(x_local,y_local+1);
    for k=1:1:2
        if(in_dir_index(j,k)==1)
            in_fic_upper_x=in_fic_for_out_ux(j,k);
            in_fic_upper_y=in_fic_for_out_uy(j,k);
            break;
        end
    end
	%compute upper out_disp_x
	x_disp=[in_fic_upper_x;d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_upper_disp_x=a*dis_x^2+b*dis_x+c;
	%compute upper out_disp_y
	y_disp=[in_fic_upper_y;d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_upper_disp_y=a*dis_x^2+b*dis_x+c;
	% upper2 situation
	j=in_fic_index(x_local,y_local+2);
    for k=1:1:2
        if(in_dir_index(j,k)==1)
            in_fic_upper2_x=in_fic_for_out_ux(j,k);
            in_fic_upper2_y=in_fic_for_out_uy(j,k);
            break;
        end
    end
	%compute upper2 out_disp_x
	x_disp=[in_fic_upper2_x;d_d(x_global,y_global+2,1);d_d(x_global+1,y_global+2,1)];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
	%compute upper2 out_disp_y            
	y_disp=[in_fic_upper2_y;d_d(x_global,y_global+2,2);d_d(x_global+1,y_global+2,2)];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	out_upper2_disp_y=a*dis_x^2+b*dis_x+c;  
	%compute y direction derivative
	out_deri_y_disp_x(size_in_fic_for_out,1)=(-3*out_disp_x(size_in_fic_for_out,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
	out_deri_y_disp_y(size_in_fic_for_out,1)=(-3*out_disp_y(size_in_fic_for_out,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
	
	for i=2:1:size_in_fic_for_out-1
		x_global=in_index_x(i);
		y_global=in_index_y(i);
		x_local=x_global;
		y_local=y_global;
		if(in_fic_for_out_ux(i,2)==0)
			str_variable=char(in_fic_for_out_ux(i,1));
			length_str=length(str_variable);
			if(str_variable(length_str-1)=='_')
				if(out_cross_boun_mark(x_local,y_local,1)==1)
                    if(str_variable(length_str)=='L')            
                        dis_y=(in_index_y(i)-1)*d_y-center_y;
                        dis_x=sqrt(r^2-dis_y^2)+center_x;
                        %the index of in_fic on outer point (x,y+1)
                        %the corresponding in_fic_point is on the left of (x,y+1)
                        %upper situation
                        j=in_fic_index(x_local,y_local+1);
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
                        x_disp=[in_fic_upper_x;d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                        %compute upper out_disp_y
                        y_disp=[in_fic_upper_y;d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                        % bottom situation
                        j=in_fic_index(x_local,y_local-1);
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
                        x_disp=[in_fic_bottom_x;d_d(x_global,y_global-1,1);d_d(x_global+1,y_global-1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                        %compute bottom out_disp_y            
                        y_disp=[in_fic_bottom_y;d_d(x_global,y_global-1,2);d_d(x_global+1,y_global-1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_bottom_disp_y=a*dis_x^2+b*dis_x+c;  
                        %compute y direction derivative
                        out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
                        out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);    
                    end
					if(str_variable(length_str)=='B')
                        dis_x=(in_index_x(i)-1)*d_x-center_x;
                        dis_y=sqrt(r^2-dis_x^2)+center_y;
                        %left situation
                        j=in_fic_index(x_local-1,y_local);
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
                        x_disp=[in_fic_left_x;d_d(x_global-1,y_global,1);d_d(x_global-1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_left_disp_x=a*dis_y^2+b*dis_y+c;
                        %compute left out_disp_y
                        y_disp=[in_fic_left_y;d_d(x_global-1,y_global,2);d_d(x_global-1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_left_disp_y=a*dis_y^2+b*dis_y+c;
                        %right situation
                        j=in_fic_index(x_local+1,y_local);
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
                        x_disp=[in_fic_right_x;d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_x=a*dis_y^2+b*dis_y+c;
                        %compute right out_disp_y
                        y_disp=[in_fic_right_y;d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_y=a*dis_y^2+b*dis_y+c;
                        out_deri_x_disp_x(i,1)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
                        out_deri_x_disp_y(i,1)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
					end
				elseif(out_cross_boun_mark(x_local,y_local,1)==4)
                    if(str_variable(length_str)=='L')
                        dis_y=(in_index_y(i)-1)*d_y-center_y;
                        dis_x=sqrt(r^2-dis_y^2)+center_x;
                        %upper points, real points
                        %compute upper out_disp_x
                        x_disp=[d_d(x_global-1,y_global+1,1);d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                        %compute upper out_disp_y
                        y_disp=[d_d(x_global-1,y_global+1,2);d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                        % upper2 situation
                        x_disp=[d_d(x_global-1,y_global+2,1);d_d(x_global,y_global+2,1);d_d(x_global+1,y_global+2,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
                        %compute upper2 out_disp_y
                        y_disp=[d_d(x_global-1,y_global+2,2);d_d(x_global,y_global+2,2);d_d(x_global+1,y_global+2,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
                        %compute y direction derivative
                        out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
                        out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y); 
                    end
                    if(str_variable(length_str)=='B')
                        dis_x=(in_index_x(i)-1)*d_x-center_x;
                        dis_y=sqrt(r^2-dis_x^2)+center_y;
                        %right2 situation
                        %compute right out_disp_x
                        x_disp=[d_d(x_global+2,y_global-1,1);d_d(x_global+2,y_global,1);d_d(x_global+2,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_x=a*dis_y^2+b*dis_y+c;
                        %compute right out_disp_y
                        y_disp=[d_d(x_global+2,y_global-1,2);d_d(x_global+2,y_global,2);d_d(x_global+2,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_y=a*dis_y^2+b*dis_y+c;
                        %right situation,real point
                        %compute right out_disp_x
                        x_disp=[d_d(x_global+1,y_global-1,1);d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_x=a*dis_y^2+b*dis_y+c;
                        %compute right out_disp_y
                        y_disp=[d_d(x_global+1,y_global-1,2);d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_y=a*dis_y^2+b*dis_y+c;
                        out_deri_x_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
                        out_deri_x_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
                    end
				end
			else
				continue;
			end
		elseif(in_fic_for_out_ux(i,2)~=0)
			if(out_cross_boun_mark(x_local,y_local,1)==2) 
            %left situation:('L')
                dis_y=(in_index_y(i)-1)*d_y-center_y;
                dis_x=sqrt(r^2-dis_y^2)+center_x;
                %upper situation
                %ux
                x_disp=[d_d(x_global-1,y_global+1,1);d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                %uy
                y_disp=[d_d(x_global-1,y_global+1,2);d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                %bottom situation
                %ux
                x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,2);d_d(x_global+1,y_global-1,1)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                %uy
                y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,2);d_d(x_global+1,y_global-1,2)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
                out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
            %bottom situatio:('B')
                dis_x=(in_index_x(i)-1)*d_x-center_x;
                dis_y=sqrt(r^2-dis_x^2)+center_y;
                %right situation
                %ux
                x_disp=[d_d(x_global+1,y_global-1,1);d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_x=a*dis_y^2+b*dis_y+c;
                %uy
                y_disp=[d_d(x_global+1,y_global-1,2);d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_right_disp_y=a*dis_y^2+b*dis_y+c;
                %left situation
                %ux
                x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,1);d_d(x_global-1,y_global+1,1)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_left_disp_x=a*dis_y^2+b*dis_y+c;
                %uy
                y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,1);d_d(x_global-1,y_global+1,2)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                out_left_disp_y=a*dis_y^2+b*dis_y+c;
                out_deri_x_disp_x(i,2)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
                out_deri_x_disp_y(i,2)=(out_right_disp_y-out_left_disp_y)/(2*d_x);
			elseif(out_cross_boun_mark(x_local,y_local,1)==5)
                if(out_cross_boun_mark(x_local,y_local,2)==1)
                %left situation:('L')
                    dis_y=(in_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    %upper situation
                    %ux
                    x_disp=[d_d(x_global-1,y_global+1,1);d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                    y_disp=[d_d(x_global-1,y_global+1,2);d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                    %bottom situation
                    %ux
                    j=in_fic_index(x_local+1,y_local);
                    x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,2);in_fic_for_out_ux(j,1)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                    %uy
                    y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,2);in_fic_for_out_uy(j,1)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                    out_deri_y_disp_x(i,1)=(out_upper_disp_x-out_bottom_disp_x)/(2*d_y);
                    out_deri_y_disp_y(i,1)=(out_upper_disp_y-out_bottom_disp_y)/(2*d_y);
                %bottom situatio:('B')
                    dis_x=(in_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y;
                    %right situation
                    j_1=in_fic_index(x_local+1,y_local);
                    j_2=in_fic_index(x_local+2,y_local);
                    if(local_point(x_local+2,y_local-1)==-1)
                        x_disp=[in_fic_for_out_ux(j_1,1);d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[in_fic_for_out_uy(j_1,1);d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
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
                        x_disp=[in_right2_x;d_d(x_global+2,y_global,1);d_d(x_global+2,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[in_right2_y;d_d(x_global+2,y_global,2);d_d(x_global+2,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_y=a*dis_y^2+b*dis_y+c;
                        out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
                        out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
                        elseif(local_point(x_local+2,y_local-1)==1)
                            for k=1:1:2
                                if(in_dir_index(j_1,k)==2)
                                    in_right_x=in_fic_for_out_ux(j_1,k);
                                    in_right_y=in_fic_for_out_uy(j_1,k);
                                end
                            end
                        x_disp=[in_right_x;d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[in_right_y;d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right_disp_y=a*dis_y^2+b*dis_y+c;
                        x_disp=[d_d(x_global+2,y_global-1,1);d_d(x_global+2,y_global,1);d_d(x_global+2,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[d_d(x_global+2,y_global-1,2);d_d(x_global+2,y_global,2);d_d(x_global+2,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_right2_disp_y=a*dis_y^2+b*dis_y+c;
                        out_deri_x_disp_x(i,2)=(-3*out_disp_x(i,2)+4*out_right_disp_x-out_right2_disp_x)/(2*d_x);
                        out_deri_x_disp_y(i,2)=(-3*out_disp_y(i,2)+4*out_right_disp_y-out_right2_disp_y)/(2*d_x);
                    end
                end
                if(out_cross_boun_mark(x_local,y_local,2)==2)
                %left situation:('L')
                    dis_y=(in_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    j_1=in_fic_index(x_local,y_local+1);
                    j_2=in_fic_index(x_local,y_local+2);
                    if(local_point(x_local-1,y_local+2)==-1)
                            %ux
                        x_disp=[in_fic_for_out_ux(j_1,1);d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                            %uy
                        y_disp=[in_fic_for_out_uy(j_1,1);d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
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
                        x_disp=[out_upper2_x;d_d(x_global,y_global+2,1);d_d(x_global+1,y_global+2,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
                            %uy
                        y_disp=[out_upper2_y;d_d(x_global,y_global+2,2);d_d(x_global+1,y_global+2,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
                        out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
                        out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
                    elseif(local_point(x_local-1,y_local+2)==1)
                        for k=1:1:2
                            if(in_dir_index(j_1,k)==1)
                                out_upper_x=in_fic_for_out_ux(j_1,k);
                                out_upper_y=in_fic_for_out_uy(j_1,k);
                            end
                        end
                            %ux
                        x_disp=[out_upper_x;d_d(x_global,y_global+1,1);d_d(x_global+1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_x=a*dis_x^2+b*dis_x+c;
                            %uy
                        y_disp=[out_upper_y;d_d(x_global,y_global+1,2);d_d(x_global+1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper_disp_y=a*dis_x^2+b*dis_x+c;
                            %ux
                        x_disp=[d_d(x_global-1,y_global+2,1);d_d(x_global,y_global+2,1);d_d(x_global+1,y_global+2,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_x=a*dis_x^2+b*dis_x+c;
                            %uy
                        y_disp=[d_d(x_global-1,y_global+2,2);d_d(x_global,y_global+2,2);d_d(x_global+1,y_global+2,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        out_upper2_disp_y=a*dis_x^2+b*dis_x+c;
                        out_deri_y_disp_x(i,1)=(-3*out_disp_x(i,1)+4*out_upper_disp_x-out_upper2_disp_x)/(2*d_y);
                        out_deri_y_disp_y(i,1)=(-3*out_disp_y(i,1)+4*out_upper_disp_y-out_upper2_disp_y)/(2*d_y);
                    end
                %bottom situatio:('B')
                    dis_x=(in_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y;
                    %right situation
                    %ux
                    x_disp=[d_d(x_global+1,y_global-1,1);d_d(x_global+1,y_global,1);d_d(x_global+1,y_global+1,1)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_right_disp_x=a*dis_y^2+b*dis_y+c;
                    %uy
                    y_disp=[d_d(x_global+1,y_global-1,2);d_d(x_global+1,y_global,2);d_d(x_global+1,y_global+1,2)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_right_disp_y=a*dis_y^2+b*dis_y+c;
                %left situation
                    %ux
                    j=in_fic_index(x_local,y_local+1);
                    x_disp=[in_fic_for_out_ux(i,3);in_fic_for_out_ux(i,1);in_fic_for_out_ux(j,1)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_left_disp_x=a*dis_y^2+b*dis_y+c;
                    %uy
                    y_disp=[in_fic_for_out_uy(i,3);in_fic_for_out_uy(i,1);in_fic_for_out_uy(j,1)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    out_left_disp_y=a*dis_y^2+b*dis_y+c;
                    out_deri_x_disp_x(i,2)=(out_right_disp_x-out_left_disp_x)/(2*d_x);
                    out_deri_x_disp_y(i,2)=(out_right_disp_y-out_left_disp_y)/(2*d_x);  
                end
			end
		end
	end
	
	% boundary points
	x_global=out_index_x(size_out_fic_for_in);
	y_global=out_index_y(size_out_fic_for_in);
	x_local=x_global;
	y_local=y_global;
	dis_y=(out_index_y(size_out_fic_for_in)-1)*d_y-center_y;
	dis_x=sqrt(r^2-dis_y^2)+center_x;
	%upper situation
	j=out_fic_index(x_local,y_local+1);
    for k=1:1:2
        if(out_dir_index(j,k)==1)
            out_fic_upper_x=out_fic_for_in_ux(j,k);
            out_fic_upper_y=out_fic_for_in_uy(j,k);
            break;
        end
    end
	%upper ux
	x_disp=[d_d(x_global-1,y_global+1,1);d_d(x_global,y_global+1,1);out_fic_upper_x];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_upper_disp_x=a*dis_x^2+b*dis_x+c;
	%upper uy            
	y_disp=[d_d(x_global-1,y_global+1,2);d_d(x_global,y_global+1,2);out_fic_upper_y];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_upper_disp_y=a*dis_x^2+b*dis_x+c;
	%upper2 situation
	j=out_fic_index(x_local,y_local+2);
    for k=1:1:2
        if(out_dir_index(j,k)==1)
            out_fic_upper2_x=out_fic_for_in_ux(j,k);
            out_fic_upper2_y=out_fic_for_in_uy(j,k);
            break;
        end
    end
	%upper2 ux
	x_disp=[d_d(x_global-1,y_global+2,1);d_d(x_global,y_global+2,1);out_fic_upper2_x];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_upper2_disp_x=a*dis_x^2+b*dis_x+c;
	%upper2 uy            
	y_disp=[d_d(x_global-1,y_global+2,2);d_d(x_global,y_global+2,2);out_fic_upper2_y];
	syms a b c;
	eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
	eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
	eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_upper2_disp_y=a*dis_x^2+b*dis_x+c;
	in_deri_y_disp_x(size_out_fic_for_in,1)=(-3*in_disp_x(size_out_fic_for_in,1)+4*in_upper_disp_x-in_upper2_disp_x)/(2*d_y);
	in_deri_y_disp_y(size_out_fic_for_in,1)=(-3*in_disp_y(size_out_fic_for_in,1)+4*in_upper_disp_y-in_upper2_disp_y)/(2*d_y);
	
	%boundary points
	x_global=out_index_x(1);
	y_global=out_index_y(1);
	x_local=x_global;
	y_local=y_global;
	dis_x=(out_index_x(1)-1)*d_x-center_x;
	dis_y=sqrt(r^2-dis_x^2)+center_y; 
	%right2 situation
	j=out_fic_index(x_local+2,y_local);
    for k=1:1:2
        if(out_dir_index(j,k)==2)
            out_fic_right2_ux=out_fic_for_in_ux(j,k);
            out_fic_right2_uy=out_fic_for_in_uy(j,k);
            break;
        end
    end
	%right2 ux
	x_disp=[d_d(x_global+2,y_global-1,1);d_d(x_global+2,y_global,1);out_fic_right2_ux];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_right2_disp_x=a*dis_y^2+b*dis_y+c;
	%right2 uy
	y_disp=[d_d(x_global+2,y_global-1,2);d_d(x_global+2,y_global,2);out_fic_right2_uy];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_right2_disp_y=a*dis_y^2+b*dis_y+c;
	%right situation
	j=out_fic_index(x_local+1,y_local);
    for k=1:1:2
        if(out_dir_index(j,k)==2)
            out_fic_right_ux=out_fic_for_in_ux(j,k);
            out_fic_right_uy=out_fic_for_in_uy(j,k);
            break;
        end
    end
	%right ux
	x_disp=[d_d(x_global+1,y_global-1,1);d_d(x_global+1,y_global,1);out_fic_right_ux];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_right_disp_x=a*dis_y^2+b*dis_y+c;
	%right uy
	y_disp=[d_d(x_global+1,y_global-1,2);d_d(x_global+1,y_global,2);out_fic_right_uy];
	syms a b c;
	eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
	eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
	eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
	[a,b,c]=solve(eqn,[a,b,c]);
	in_right_disp_y=a*dis_y^2+b*dis_y+c; 
	in_deri_x_disp_x(1,1)=(-3*in_disp_x(1,1)+4*in_right_disp_x-in_right2_disp_x)/(2*d_x);
	in_deri_x_disp_y(1,1)=(-3*in_disp_y(1,1)+4*in_right_disp_y-in_right2_disp_y)/(2*d_x);
	
    for i=2:1:size_out_fic_for_in-1
        x_global=out_index_x(i);
        y_global=out_index_y(i);
        x_local=x_global;
        y_local=y_global;
        if(out_fic_for_in_ux(i,2)==0)
            str_variable=char(out_fic_for_in_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1)=='_')
                if(in_cross_boun_mark(x_local,y_local,1)==1)
                    if(str_variable(length_str)=='R')
                    dis_y=(out_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    %upper situation
                    j=out_fic_index(x_local,y_local+1);
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
                    x_disp=[d_d(x_global-1,y_global+1,1);d_d(x_global,y_global+1,1);out_fic_upper_x];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_upper_disp_x=a*dis_x^2+b*dis_x+c;
                    %upper uy            
                    y_disp=[d_d(x_global-1,y_global+1,2);d_d(x_global,y_global+1,2);out_fic_upper_y];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_upper_disp_y=a*dis_x^2+b*dis_x+c;
                    %bottom situation
                    j=out_fic_index(x_local,y_local-1);
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
                    x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);out_fic_bottom_x];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_bottom_disp_x=a*dis_x^2+b*dis_x+c;   
                    %bottom uy
                    y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);out_fic_bottom_y];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_bottom_disp_y=a*dis_x^2+b*dis_x+c;  
                    in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
                    in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
                    end
                    if(str_variable(length_str)=='T')
                    dis_x=(out_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y; 
                    %left situation
                    j=out_fic_index(x_local-1,y_local);
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
                    x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);out_fic_left_ux];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_left_disp_x=a*dis_y^2+b*dis_y+c;
                    %left uy
                    y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);out_fic_left_uy];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_left_disp_y=a*dis_y^2+b*dis_y+c;
                    %right situation
                    j=out_fic_index(x_local+1,y_local);
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
                    x_disp=[d_d(x_global+1,y_global-1,1);d_d(x_global+1,y_global,1);out_fic_right_ux];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_right_disp_x=a*dis_y^2+b*dis_y+c;
                    %right uy
                    y_disp=[d_d(x_global+1,y_global-1,2);d_d(x_global+1,y_global,2);out_fic_right_uy];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_right_disp_y=a*dis_y^2+b*dis_y+c; 
                    in_deri_x_disp_x(i,1)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
                    in_deri_x_disp_y(i,1)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
                    end
                elseif(in_cross_boun_mark(x_local,y_local,1)==5)
                    if(in_cross_boun_mark(x_local,y_local,2)==1)
                        dis_x=(out_index_x(i)-1)*d_x-center_x;
                        dis_y=sqrt(r^2-dis_x^2)+center_y;
                        %left2 situation
                        x_disp=[d_d(x_global-2,y_global-1,1);d_d(x_global-2,y_global,1);d_d(x_global-2,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_x=a*dis_y^2+b*dis_y+c;
                        %uy
                        y_disp=[d_d(x_global-2,y_global-1,2);d_d(x_global-2,y_global,2);d_d(x_global-2,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_y=a*dis_y^2+b*dis_y+c;
                        %left situation
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);d_d(x_global-1,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left_disp_x=a*dis_y^2+b*dis_y+c;
                        %uy
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);d_d(x_global-1,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left_disp_y=a*dis_y^2+b*dis_y+c;
                        in_deri_x_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
                        in_deri_x_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
                    end
                    if(in_cross_boun_mark(x_local,y_local,2)==2)
                        dis_y=(out_index_y(i)-1)*d_y-center_y;
                        dis_x=sqrt(r^2-dis_y^2)+center_x;
                        %bottom2 situation
                        x_disp=[d_d(x_global-1,y_global-2,1);d_d(x_global,y_global-2,1);d_d(x_global+1,y_global-2,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
                        %uy
                        y_disp=[d_d(x_global-1,y_global-2,2);d_d(x_global,y_global-2,2);d_d(x_global+1,y_global-2,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;  
                        %bottom situation
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);d_d(x_global+1,y_global-1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                        %uy
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);d_d(x_global+1,y_global-1,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                        in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
                        in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
                    end
                end
            else
                continue;
            end
        elseif(out_fic_for_in_ux(i,2)~=0)
            if(in_cross_boun_mark(x_local,y_local,1)==3)
            %right situation:('R')
                dis_y=(out_index_y(i)-1)*d_y-center_y;
                dis_x=sqrt(r^2-dis_y^2)+center_x;
                %bottom situation
                x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);d_d(x_global+1,y_global-1,1)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);d_d(x_global+1,y_global-1,2)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                %top situation
                x_disp=[d_d(x_global-1,y_global+1,1);out_fic_for_in_ux(i,2);out_fic_for_in_ux(i,3)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_upper_disp_x=a*dis_x^2+b*dis_x+c;
                y_disp=[d_d(x_global-1,y_global+1,2);out_fic_for_in_uy(i,2);out_fic_for_in_uy(i,3)];
                syms a b c;
                eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_upper_disp_y=a*dis_x^2+b*dis_x+c;
                in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
                in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
            %upper situation:('T')
                dis_x=(out_index_x(i)-1)*d_x-center_x;
                dis_y=sqrt(r^2-dis_x^2)+center_y;
                %left situation
                x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);d_d(x_global-1,y_global+1,1)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_x=a*dis_y^2+b*dis_y+c;  
                y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);d_d(x_global-1,y_global+1,2)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_left_disp_y=a*dis_y^2+b*dis_y+c;  
                %right situation
                x_disp=[d_d(x_global+1,y_global-1,1);out_fic_for_in_ux(i,1);out_fic_for_in_ux(i,3)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_right_disp_x=a*dis_y^2+b*dis_y+c;  
                y_disp=[d_d(x_global+1,y_global-1,2);out_fic_for_in_uy(i,1);out_fic_for_in_uy(i,3)];
                syms a b c;
                eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                [a,b,c]=solve(eqn,[a,b,c]);
                in_right_disp_y=a*dis_y^2+b*dis_y+c;
                in_deri_x_disp_x(i,2)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
                in_deri_x_disp_y(i,2)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
            elseif(in_cross_boun_mark(x_local,y_local,1)==4)
                if(in_cross_boun_mark(x_local,y_local,2)==1)
                %right situation:('R')
                    dis_y=(out_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    %bottom
                    x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);d_d(x_global+1,y_global-1,1)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                    y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);d_d(x_global+1,y_global-1,2)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_bottom_disp_y=a*dis_x^2+b*dis_x+c;   
                    %upper
                    j=out_fic_index(x_local-1,y_local);
                    x_disp=[out_fic_for_in_ux(j,1);out_fic_for_in_ux(i,2);out_fic_for_in_ux(i,3)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_upper_disp_x=a*dis_x^2+b*dis_x+c;
                    y_disp=[out_fic_for_in_uy(j,1);out_fic_for_in_uy(i,2);out_fic_for_in_uy(i,3)];
                    syms a b c;
                    eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                    eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                    eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_upper_disp_y=a*dis_x^2+b*dis_x+c;
                    in_deri_y_disp_x(i,1)=(in_upper_disp_x-in_bottom_disp_x)/(2*d_y);
                    in_deri_y_disp_y(i,1)=(in_upper_disp_y-in_bottom_disp_y)/(2*d_y);
                %upper situation:('T')
                    dis_x=(out_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y; 
                    %left
                    j_1=out_fic_index(x_local-1,y_local);
                    j_2=out_fic_index(x_local-2,y_local);
                    if(local_point(x_local-2,y_local+1)==1)
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);out_fic_for_in_ux(j_1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);out_fic_for_in_uy(j_1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
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
                        x_disp=[d_d(x_global-2,y_global-1,1);d_d(x_global-2,y_global,1);out_left2_x];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[d_d(x_global-2,y_global-1,2);d_d(x_global-2,y_global,2);out_left2_y];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_y=a*dis_y^2+b*dis_y+c;
                        in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
                        in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
                    elseif(local_point(x_local-2,y_local+1)==-1)
                        for k=1:1:2
                            if(out_dir_index(j_1,k)==2)
                                out_left_x=out_fic_for_in_ux(j_1,k);
                                out_left_y=out_fic_for_in_uy(j_1,k);
                            end
                        end
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);out_left_x];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);out_left_y];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left_disp_y=a*dis_y^2+b*dis_y+c;
                        x_disp=[d_d(x_global-2,y_global-1,1);d_d(x_global-2,y_global,1);d_d(x_global-2,y_global+1,1)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_x=a*dis_y^2+b*dis_y+c;
                        y_disp=[d_d(x_global-2,y_global-1,2);d_d(x_global-2,y_global,2);d_d(x_global-2,y_global+1,2)];
                        syms a b c;
                        eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                        eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                        eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_left2_disp_y=a*dis_y^2+b*dis_y+c;
                        in_deri_x_disp_x(i,2)=(3*in_disp_x(i,2)-4*in_left_disp_x+in_left2_disp_x)/(2*d_x);
                        in_deri_x_disp_y(i,2)=(3*in_disp_y(i,2)-4*in_left_disp_y+in_left2_disp_y)/(2*d_x);
                    end
                end
                if(in_cross_boun_mark(x_local,y_local,2)==2)
                %right situation:('R')
                    dis_y=(out_index_y(i)-1)*d_y-center_y;
                    dis_x=sqrt(r^2-dis_y^2)+center_x;
                    %bottom
                    j_1=out_fic_index(x_local,y_local-1);
                    j_2=out_fic_index(x_local,y_local-2);
                    if(local_point(x_local+1,y_local-2)==1)
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);out_fic_for_in_ux(j_1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);out_fic_for_in_uy(j_1,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
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
                        x_disp=[d_d(x_global-1,y_global-2,1);d_d(x_global,y_global-2,1);out_bottom2_x];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
                        y_disp=[d_d(x_global-1,y_global-2,2);d_d(x_global,y_global-2,2);out_bottom2_y];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
                        in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
                        in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
                    elseif(local_point(x_local+1,y_local-2)==-1)
                        for k=1:1:2
                            if(out_dir_index(j_1,k)==1)
                                out_bottom_x=out_fic_for_in_ux(j_1,k);
                                out_bottom_y=out_fic_for_in_uy(j_1,k);
                            end
                        end
                        x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global,y_global-1,1);out_bottom_x];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom_disp_x=a*dis_x^2+b*dis_x+c;
                        y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global,y_global-1,2);out_bottom_y];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom_disp_y=a*dis_x^2+b*dis_x+c;
                        x_disp=[d_d(x_global-1,y_global-2,1);d_d(x_global,y_global-2,1);d_d(x_global+1,y_global-2,1)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==x_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==x_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==x_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_x=a*dis_x^2+b*dis_x+c;
                        y_disp=[d_d(x_global-1,y_global-2,2);d_d(x_global,y_global-2,2);d_d(x_global+1,y_global-2,2)];
                        syms a b c;
                        eqn(1)=a*((x_global-2)*d_x)^2+b*((x_global-2)*d_x)+c==y_disp(1);
                        eqn(2)=a*((x_global-1)*d_x)^2+b*((x_global-1)*d_x)+c==y_disp(2);
                        eqn(3)=a*((x_global)*d_x)^2+b*((x_global)*d_x)+c==y_disp(3);
                        [a,b,c]=solve(eqn,[a,b,c]);
                        in_bottom2_disp_y=a*dis_x^2+b*dis_x+c;
                        in_deri_y_disp_x(i,1)=(3*in_disp_x(i,1)-4*in_bottom_disp_x+in_bottom2_disp_x)/(2*d_y);
                        in_deri_y_disp_y(i,1)=(3*in_disp_y(i,1)-4*in_bottom_disp_y+in_bottom2_disp_y)/(2*d_y);
                    end
                %upper situation:('T') 
                    dis_x=(out_index_x(i)-1)*d_x-center_x;
                    dis_y=sqrt(r^2-dis_x^2)+center_y; 
                    %left
                    x_disp=[d_d(x_global-1,y_global-1,1);d_d(x_global-1,y_global,1);d_d(x_global-1,y_global+1,1)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_left_disp_x=a*dis_y^2+b*dis_y+c;
                    y_disp=[d_d(x_global-1,y_global-1,2);d_d(x_global-1,y_global,2);d_d(x_global-1,y_global+1,2)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_left_disp_y=a*dis_y^2+b*dis_y+c;
                    %right
                    j=out_fic_index(x_local,y_local-1);
                    x_disp=[out_fic_for_in_ux(j,1);out_fic_for_in_ux(i,1);out_fic_for_in_ux(i,3)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==x_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==x_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==x_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_right_disp_x=a*dis_y^2+b*dis_y+c;
                    y_disp=[out_fic_for_in_uy(j,1);out_fic_for_in_uy(i,1);out_fic_for_in_uy(i,3)];
                    syms a b c;
                    eqn(1)=a*((y_global-2)*d_y)^2+b*((y_global-2)*d_y)+c==y_disp(1);
                    eqn(2)=a*((y_global-1)*d_y)^2+b*((y_global-1)*d_y)+c==y_disp(2);
                    eqn(3)=a*((y_global)*d_y)^2+b*((y_global)*d_y)+c==y_disp(3);
                    [a,b,c]=solve(eqn,[a,b,c]);
                    in_right_disp_y=a*dis_y^2+b*dis_y+c;
                    in_deri_x_disp_x(i,2)=(in_right_disp_x-in_left_disp_x)/(2*d_x);
                    in_deri_x_disp_y(i,2)=(in_right_disp_y-in_left_disp_y)/(2*d_x);
                end
            end
        end
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
    
	in_fic_vari_sqn=zeros(2*length(out_deri_x_disp_x),1);
	in_fic_variables=sym(zeros(2*length(out_deri_x_disp_x)));
	out_fic_variables=sym(zeros(2*length(out_deri_x_disp_x)));
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
	
	out_sur_nor_dir_x=out_sur_nor_dir_x/r;
	out_sur_nor_dir_y=out_sur_nor_dir_y/r;
	%% discontinuity boundary condition setting
	M_steel=2*miu_steel*(1-v_steel)/(1-2*v_steel);
	M_aluminum=2*miu_aluminum*(1-v_aluminum)/(1-2*v_aluminum);

	syms sin_theta cos_theta;
    boun_traction_1=[M_aluminum*cos_theta,-M_steel*cos_theta,miu_aluminum*sin_theta,-miu_steel*sin_theta,...
        miu_aluminum*sin_theta,-miu_steel*sin_theta,M_aluminum*v_aluminum*cos_theta/(1-v_aluminum),-lambda_steel*v_steel*cos_theta/(1-v_steel)];
    boun_traction_2=[M_aluminum*v_aluminum*sin_theta/(1-v_aluminum),-lambda_steel*v_steel*sin_theta/(1-v_steel),miu_aluminum*cos_theta,-miu_steel*cos_theta,...
        miu_aluminum*cos_theta,-miu_steel*cos_theta,M_aluminum*sin_theta,-M_steel*sin_theta];

    
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
            num_boun_trac_1=subs(boun_traction_1,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            num_boun_trac_2=subs(boun_traction_2,[sin_theta,cos_theta],[out_sur_nor_dir_y(i,1),out_sur_nor_dir_x(i,1)]);
            if(in_fic_for_out_ux(j,2)==0)
                C=[in_deri_x_disp_x(i,1);out_deri_x_disp_x(j,1);in_deri_y_disp_x(i,1);out_deri_y_disp_x(j,1);...
                    in_deri_x_disp_y(i,1);out_deri_x_disp_y(j,1);in_deri_y_disp_y(i,1);out_deri_y_disp_y(j,1)];
                disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,1)==0;
                disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,1)==0;
                disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                index_eqn=index_eqn+4;
            else
                for m=1:1:2
                    if(out_dir_index(i,1)==1)
                        if(in_dir_index(j,m)==1)
                        C=[in_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);in_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                            in_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);in_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
                        disconti_boun_eqn(index_eqn)=in_disp_x(i,1)-out_disp_x(j,m)==0;
                        disconti_boun_eqn(index_eqn+1)=in_disp_y(i,1)-out_disp_y(j,m)==0;
                        disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                        disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                        index_eqn=index_eqn+4;
                        end
                    elseif(out_dir_index(i,1)==2)
                        if(in_dir_index(j,m)==2)
                        C=[in_deri_x_disp_x(i,1);out_deri_x_disp_x(j,m);in_deri_y_disp_x(i,1);out_deri_y_disp_x(j,m);...
                            in_deri_x_disp_y(i,1);out_deri_x_disp_y(j,m);in_deri_y_disp_y(i,1);out_deri_y_disp_y(j,m)];
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
                    C=[in_deri_x_disp_x(i,k);out_deri_x_disp_x(j,1);in_deri_y_disp_x(i,k);out_deri_y_disp_x(j,1);...
                        in_deri_x_disp_y(i,k);out_deri_x_disp_y(j,1);in_deri_y_disp_y(i,k);out_deri_y_disp_y(j,1)];
                    disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,1)==0;
                    disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,1)==0;
                    disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                    disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                    index_eqn=index_eqn+4;
                else
                    for m=1:1:2
                        if(out_dir_index(i,k)==1)
                            if(in_dir_index(j,m)==1)
                            C=[in_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);in_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                in_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);in_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
                            disconti_boun_eqn(index_eqn)=in_disp_x(i,k)-out_disp_x(j,m)==0;
                            disconti_boun_eqn(index_eqn+1)=in_disp_y(i,k)-out_disp_y(j,m)==0;
                            disconti_boun_eqn(index_eqn+2)=num_boun_trac_1*C==0;
                            disconti_boun_eqn(index_eqn+3)=num_boun_trac_2*C==0;
                            index_eqn=index_eqn+4;
                            end
                        elseif(out_dir_index(i,k)==2)
                            if(in_dir_index(j,m)==2)
                            C=[in_deri_x_disp_x(i,k);out_deri_x_disp_x(j,m);in_deri_y_disp_x(i,k);out_deri_y_disp_x(j,m);...
                                in_deri_x_disp_y(i,k);out_deri_x_disp_y(j,m);in_deri_y_disp_y(i,k);out_deri_y_disp_y(j,m)];
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
	
	for i=1:1:length(disconti_boun_eqn)
		if(disconti_boun_eqn(i)==0)
			disconti_boun_eqn=disconti_boun_eqn(1:1:i-1);
			break;
		end
	end
	disconti_boun_eqn=vpa(disconti_boun_eqn);
	%% solve disconti_boun_eqns
	fic_variables=[in_fic_variables,out_fic_variables];
	fic_variables=transpose(fic_variables);
	solutions=solve(disconti_boun_eqn,fic_variables);
	fields=fieldnames(solutions);
	in_fic_vari_x=sym(zeros(length(in_fic_for_out_ux),3));
	in_fic_vari_y=sym(zeros(length(in_fic_for_out_ux),3));
	out_fic_vari_x=sym(zeros(length(out_fic_for_in_ux),3));
	out_fic_vari_y=sym(zeros(length(out_fic_for_in_ux),3));
	
    for i=1:1:length(in_fic_vari_sqn)
        j=in_fic_vari_sqn(i);
        if(i==1)
            if(j~=in_fic_vari_sqn(i+1))
                in_fic_vari_x(j,1)=solutions.(fields{2*i-1});
                in_fic_vari_y(j,1)=solutions.(fields{2*i});
            elseif(j==in_fic_vari_sqn(i+1))
                in_fic_vari_x(j,1)=solutions.(fields{2*i-1});
                in_fic_vari_y(j,1)=solutions.(fields{2*i});
                in_fic_vari_x(j,2)=solutions.(fields{2*i+1});
                in_fic_vari_y(j,2)=solutions.(fields{2*i+2});
            end
        elseif(i<length(in_fic_vari_sqn))
            if(j==in_fic_vari_sqn(i-1))
                continue;
            else
                if(j~=in_fic_vari_sqn(i+1))
                    in_fic_vari_x(j,1)=solutions.(fields{2*i-1});
                    in_fic_vari_y(j,1)=solutions.(fields{2*i});
                elseif(j==in_fic_vari_sqn(i+1))
                    in_fic_vari_x(j,1)=solutions.(fields{2*i-1});
                    in_fic_vari_y(j,1)=solutions.(fields{2*i});
                    in_fic_vari_x(j,2)=solutions.(fields{2*i+1});
                    in_fic_vari_y(j,2)=solutions.(fields{2*i+2});
                end
            end
        else
            if(j==in_fic_vari_sqn(i-1))
                continue;
            else
                in_fic_vari_x(j,1)=solutions.(fields{2*i-1});
                in_fic_vari_y(j,1)=solutions.(fields{2*i});
            end
        end
    end
    for i=1:1:length(out_fic_vari_sqn)
        j=out_fic_vari_sqn(i);
        if(i==1)
            if(j~=out_fic_vari_sqn(i+1))
                out_fic_vari_x(j,1)=solutions.(fields{2*i+length_variables-1});
                out_fic_vari_y(j,1)=solutions.(fields{2*i+length_variables});
            elseif(j==out_fic_vari_sqn(i+1))
                out_fic_vari_x(j,1)=solutions.(fields{2*i+length_variables-1});
                out_fic_vari_y(j,1)=solutions.(fields{2*i+length_variables});
                out_fic_vari_x(j,2)=solutions.(fields{2*i+length_variables+1});
                out_fic_vari_y(j,2)=solutions.(fields{2*i+length_variables+2});
            end
        elseif(i<length(out_fic_vari_sqn))
            if(j==out_fic_vari_sqn(i-1))
                continue;
            else
                if(j~=out_fic_vari_sqn(i+1))
                    out_fic_vari_x(j,1)=solutions.(fields{2*i+length_variables-1});
                    out_fic_vari_y(j,1)=solutions.(fields{2*i+length_variables});
                elseif(j==out_fic_vari_sqn(i+1))
                    out_fic_vari_x(j,1)=solutions.(fields{2*i+length_variables-1});
                    out_fic_vari_y(j,1)=solutions.(fields{2*i+length_variables});
                    out_fic_vari_x(j,2)=solutions.(fields{2*i+length_variables+1});
                    out_fic_vari_y(j,2)=solutions.(fields{2*i+length_variables+2});
                end
            end
        else
            if(j==out_fic_vari_sqn(i-1))
                continue;
            else
                out_fic_vari_x(j,1)=solutions.(fields{2*i+length_variables-1});
                out_fic_vari_y(j,1)=solutions.(fields{2*i+length_variables});
            end
        end
    end
	
    for i=1:1:length(in_fic_for_out_ux)
        index_x=in_index_x(i);
        index_y=in_index_y(i);
        if(in_fic_for_out_ux(i,2)==0)
            str_variable=char(in_fic_for_out_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1:1:length_str)=='LB')
                in_fic_vari_x(i,1)=d_d(index_x-1,index_y,1)+d_d(index_x,index_y-1,1)-d_d(index_x,index_y,1);
                in_fic_vari_y(i,1)=d_d(index_x-1,index_y,2)+d_d(index_x,index_y-1,2)-d_d(index_x,index_y,2);
            end
        end
        if(in_fic_for_out_ux(i,3)~=0)
            in_fic_vari_x(i,3)=in_fic_vari_x(i,1)+in_fic_vari_x(i,2)-d_d(index_x,index_y,1);
            in_fic_vari_y(i,3)=in_fic_vari_y(i,1)+in_fic_vari_y(i,2)-d_d(index_x,index_y,2);
        end
    end
    for i=1:1:length(out_fic_for_in_ux)
        index_x=out_index_x(i);
        index_y=out_index_y(i);
        if(out_fic_for_in_ux(i,2)==0)
            str_variable=char(out_fic_for_in_ux(i,1));
            length_str=length(str_variable);
            if(str_variable(length_str-1:1:length_str)=='RT')
                out_fic_vari_x(i,1)=d_d(index_x+1,index_y,1)+d_d(index_x,index_y+1,1)-d_d(index_x,index_y,1);
                out_fic_vari_y(i,1)=d_d(index_x+1,index_y,2)+d_d(index_x,index_y+1,2)-d_d(index_x,index_y,2);
            end
        end
        if(out_fic_for_in_ux(i,3)~=0)
            out_fic_vari_x(i,3)=out_fic_vari_x(i,1)+out_fic_vari_x(i,2)-d_d(index_x,index_y,1);
            out_fic_vari_y(i,3)=out_fic_vari_y(i,1)+out_fic_vari_y(i,2)-d_d(index_x,index_y,2);
        end
    end

    equilibrium=sym(zeros(max_x,max_y,2));
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
    %{
	for i=1:1:1
		for j=1:1:1
			equilibrium(i,j,1)=0;
			equilibrium(i,j,2)=0;
			d_d(i,j,1)=0;
			d_d(i,j,2)=0;
		end
	end
    for i=1:1:1
		for j=2:1:high_y-2
			equilibrium(i,j,1)=0;
			equilibrium(i,j,2)=miu_aluminum*(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                miu_aluminum*(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
			d_d(i,j,1)=0;
		end
	end
	for i=1:1:1
		for j=high_y-1:1:high_y-1
			equilibrium(i,j,1)=0;
			equilibrium(i,j,2)=miu_aluminum*(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                miu_aluminum*(out_fic_vari_x(1,1)-d_d(i,j-1,1))/(2*d_y);
			d_d(i,j,1)=0;
		end
	end
	for i=1:1:1
		for j=high_y:1:high_y
			equilibrium(i,j,1)=0;
			equilibrium(i,j,2)=miu_steel*(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                miu_steel*(d_d(i,j+1,1)-in_fic_vari_x(1,1))/(2*d_y);
			d_d(i,j,1)=0;
		end
	end
	for i=1:1:1
		for j=high_y+1:1:max_y-1
			equilibrium(i,j,1)=0;
			equilibrium(i,j,2)=miu_steel*(-3*d_d(i,j,2)+4*d_d(i+1,j,2)-d_d(i+2,j,2))/(2*d_x)+...
                miu_steel*(d_d(i,j+1,1)-d_d(i,j-1,1))/(2*d_y);
			d_d(i,j,1)=0;
		end
    end
    for i=1
        for j=max_y
            d_d(i,max_y,1)=0;
            equilibrium(i,max_y,1)=0;
            equilibrium(i,max_y,2)=lambda_steel*(-3*d_d(1,max_y,1)+4*d_d(2,max_y,1)-d_d(3,max_y,1))/(2*d_x)+...
    (2*miu_steel+lambda_steel)*(-3*d_d(i,max_y,2)+4*d_d(i,max_y-1,2)-d_d(i,max_y-2,2))/(2*d_y);
        end
    end
	for i=2:1:high_x-2
		for j=1:1:1
			equilibrium(i,j,1)=miu_aluminum*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_aluminum*(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=0;
			d_d(i,j,2)=0;
		end
	end
	for i=high_x-1:1:high_x-1
		for j=1:1:1
			equilibrium(i,j,1)=miu_aluminum*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_aluminum*(out_fic_vari_y(size_out_fic_for_in,1)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=0;
			d_d(i,j,2)=0;
		end
	end
	for i=high_x:1:high_x
		for j=1:1:1
			equilibrium(i,j,1)=miu_steel*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,j,2)-in_fic_vari_y(size_in_fic_for_out,1))/(2*d_x);
			equilibrium(i,j,2)=0;
			d_d(i,j,2)=0;
		end
	end
	for i=high_x+1:1:max_x-1
		for j=1:1:1
			equilibrium(i,j,1)=miu_steel*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=0;
			d_d(i,j,2)=0;
		end
    end
    
	d_d(max_x,1,2)=0;
	equilibrium(max_x,1,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
		lambda_steel*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
	equilibrium(max_x,1,2)=0;
	
	equilibrium(max_x,max_y,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
		lambda_steel*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
	equilibrium(max_x,max_y,2)=lambda_steel*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
        (2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
		
	for j=2:1:max_y-1
		equilibrium(max_x,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,j,1)-4*d_d(max_x-1,j,1)+d_d(max_x-2,j,1))/(2*d_x)+...
			lambda_steel*(d_d(max_x,j+1,2)-d_d(max_x,j-1,2))/(2*d_y);
		equilibrium(max_x,j,2)=miu_steel*(d_d(max_x,j+1,1)-d_d(max_x,j-1,1))/(2*d_y)+...
			miu_steel*(3*d_d(max_x,j,2)-4*d_d(max_x-1,j,2)+d_d(max_x-2,j,2))/(2*d_x);
	end
	for i=2:1:max_x-1
		equilibrium(i,max_y,2)=lambda_steel*(d_d(i+1,max_y,1)-d_d(i-1,max_y,1))/(2*d_x)+...
			(2*miu_steel+lambda_steel)*(3*d_d(i,max_y,2)-4*d_d(i,max_y-1,2)+d_d(i,max_y-2,2))/(2*d_y);
		equilibrium(i,max_y,1)=miu_steel*(3*d_d(i,max_y,1)-4*d_d(i,max_y-1,1)+d_d(i,max_y-2,1))/(2*d_y)+...
			miu_steel*(d_d(i+1,max_y,2)-d_d(i-1,max_y,2))/(2*d_x);
	end
	%}
    for i=2:1:size_in_fic_for_out-1
        x_global=in_index_x(i);
        y_global=in_index_y(i);
        x_local=x_global;
        y_local=y_global;
        if(out_cross_boun_mark(x_local,y_local,1)==1)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local,y_local+1);
                if(in_fic_vari_x(j_1,2)==0)
                    in_fic_upper_x=in_fic_vari_x(j_1,1);
                    in_fic_upper_y=in_fic_vari_y(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==1)
                            in_fic_upper_x=in_fic_vari_x(j_1,k);
                            in_fic_upper_y=in_fic_vari_y(j_1,k);
                        end
                    end
                end
                j_2=in_fic_index(x_local,y_local-1);
                if(in_fic_for_out_ux(j_2,2)==0)
                    in_fic_bottom_x=in_fic_vari_x(j_2,1);
                    in_fic_bottom_y=in_fic_vari_y(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==1)
                            in_fic_bottom_x=in_fic_vari_x(j_2,k);
                            in_fic_bottom_y=in_fic_vari_y(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_x(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-in_fic_upper_y-d_d(x_global+1,y_global-1,2)+in_fic_bottom_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_vari_y(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-in_fic_upper_x-d_d(x_global+1,y_global-1,1)+in_fic_bottom_x)/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local-1,y_local);
                if(in_fic_vari_x(j_1,2)==0)
                    in_fic_left_x=in_fic_vari_x(j_1,1);
                    in_fic_left_y=in_fic_vari_y(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==2)
                            in_fic_left_x=in_fic_vari_x(j_1,k);
                            in_fic_left_y=in_fic_vari_y(j_1,k);
                        end
                    end
                end
                j_2=in_fic_index(x_local+1,y_local);
                if(in_fic_for_out_ux(j_2,2)==0)
                    in_fic_right_x=in_fic_vari_x(j_2,1);
                    in_fic_right_y=in_fic_vari_y(j_2,1);
                elseif(in_fic_for_out_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_2,k)==2)
                            in_fic_right_x=in_fic_vari_x(j_2,k);
                            in_fic_right_y=in_fic_vari_y(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_vari_x(i,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-in_fic_right_y+in_fic_left_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_y(i,1)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-in_fic_right_x+in_fic_left_x)/(4*d_x*d_y));
            end
        elseif(out_cross_boun_mark(x_local,y_local,1)==2)
            equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (in_fic_vari_x(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_steel)*(in_fic_vari_x(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_vari_y(i,3))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (in_fic_vari_y(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_steel)*(in_fic_vari_y(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_vari_x(i,3))/(4*d_x*d_y));
        elseif(out_cross_boun_mark(x_local,y_local,1)==3)
            equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_vari_y(i,1))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_vari_x(i,1))/(4*d_x*d_y));
        elseif(out_cross_boun_mark(x_local,y_local,1)==4)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local-1,y_local);
                if(in_fic_vari_x(j_1,2)==0)
                    in_fic_left_x=in_fic_vari_x(j_1,1);
                    in_fic_left_y=in_fic_vari_y(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==2)
                            in_fic_left_x=in_fic_vari_x(j_1,k);
                            in_fic_left_y=in_fic_vari_y(j_1,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_vari_x(i,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_left_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_y(i,1)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_left_x)/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local,y_local-1);
                if(in_fic_vari_x(j_1,2)==0)
                    in_fic_bottom_x=in_fic_vari_x(j_1,1);
                    in_fic_bottom_y=in_fic_vari_y(j_1,1);
                elseif(in_fic_for_out_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(in_dir_index(j_1,k)==1)
                            in_fic_bottom_x=in_fic_vari_x(j_1,k);
                            in_fic_bottom_y=in_fic_vari_y(j_1,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_x(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+in_fic_bottom_y)/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_vari_y(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+in_fic_bottom_x)/(4*d_x*d_y));
            end
        elseif(out_cross_boun_mark(x_local,y_local,1)==5)
            if(out_cross_boun_mark(x_local,y_local,2)==1)
                j_1=in_fic_index(x_local+1,y_local);
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_x(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_vari_x(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-d_d(x_global-1,y_global+1,2)-in_fic_vari_y(j_1,1)+in_fic_vari_y(i,3))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_y(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_vari_y(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-d_d(x_global-1,y_global+1,1)-in_fic_vari_x(j_1,1)+in_fic_vari_x(i,3))/(4*d_x*d_y));
            elseif(out_cross_boun_mark(x_local,y_local,2)==2)
                j_1=in_fic_index(x_local,y_local+1);
                equilibrium(x_global,y_global,1)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_x(i,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_steel)*(in_fic_vari_x(i,2)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (d_d(x_global+1,y_global+1,2)-in_fic_vari_y(j_1,1)-d_d(x_global+1,y_global-1,2)+in_fic_vari_y(i,3))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_steel+lambda_steel)*(2*(1-v_steel)*...
                    (in_fic_vari_y(i,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_steel)*(in_fic_vari_y(i,1)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (d_d(x_global+1,y_global+1,1)-in_fic_vari_x(j_1,1)-d_d(x_global+1,y_global-1,1)+in_fic_vari_x(i,3))/(4*d_x*d_y));
            end

        end
    end
    for i=2:1:size_out_fic_for_in-1
        x_global=out_index_x(i);
        y_global=out_index_y(i);
        x_local=x_global;
        y_local=y_global;
        if(in_cross_boun_mark(x_local,y_local,1)==1)
           if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local,y_local+1);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_upper_x=out_fic_vari_x(j_1,1);
                    out_fic_upper_y=out_fic_vari_y(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==1)
                            out_fic_upper_x=out_fic_vari_x(j_1,k);
                            out_fic_upper_y=out_fic_vari_y(j_1,k);
                        end
                    end
                end
                j_2=out_fic_index(x_local,y_local-1);
                if(out_fic_for_in_ux(j_2,2)==0)
                    out_fic_bottom_x=out_fic_vari_x(j_2,1);
                    out_fic_bottom_y=out_fic_vari_y(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==1)
                            out_fic_bottom_x=out_fic_vari_x(j_2,k);
                            out_fic_bottom_y=out_fic_vari_y(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_upper_y-d_d(x_global-1,y_global+1,2)-out_fic_bottom_y+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_upper_x-d_d(x_global-1,y_global+1,1)-out_fic_bottom_x+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2) 
                j_1=out_fic_index(x_local-1,y_local);
                if(out_fic_for_in_ux(j_1,2)==0)
                    out_fic_left_x=out_fic_vari_x(j_1,1);
                    out_fic_left_y=out_fic_vari_y(j_1,1);
                elseif(out_fic_for_in_ux(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==2)
                            out_fic_left_x=out_fic_vari_x(j_1,k);
                            out_fic_left_y=out_fic_vari_y(j_1,k);
                        end
                    end
                end
                j_2=out_fic_index(x_local+1,y_local);
                if(out_fic_for_in_ux(j_2,2)==0)
                    out_fic_right_x=out_fic_vari_x(j_2,1);
                    out_fic_right_y=out_fic_vari_y(j_2,1);
                elseif(out_fic_for_in_ux(j_2,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_2,k)==2)
                            out_fic_right_x=out_fic_vari_x(j_2,k);
                            out_fic_right_y=out_fic_vari_y(j_2,k);
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_right_y-out_fic_left_y-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_right_x-out_fic_left_x-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        elseif(in_cross_boun_mark(x_local,y_local,1)==2)
            equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (out_fic_vari_y(i,1)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (out_fic_vari_x(i,1)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
        elseif(in_cross_boun_mark(x_local,y_local,1)==3)
            equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global-1,y_global,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_vari_x(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                (out_fic_vari_y(i,3)-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
            equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                (d_d(x_global,y_global-1,2)+out_fic_vari_y(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                (out_fic_vari_x(i,3)-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
        elseif(in_cross_boun_mark(x_local,y_local,1)==4)
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local-1,y_local);	
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_vari_x(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_vari_y(i,3)-out_fic_vari_y(j_1,1)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_vari_y(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_vari_x(i,3)-out_fic_vari_x(j_1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                j_1=out_fic_index(x_local,y_local-1);
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_vari_x(i,2)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_vari_y(i,3)-d_d(x_global-1,y_global+1,2)-out_fic_vari_y(j_1,1)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_vari_y(i,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_vari_x(i,3)-d_d(x_global-1,y_global+1,1)-out_fic_vari_x(j_1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        elseif(in_cross_boun_mark(x_local,y_local,1)==5)
            if(in_cross_boun_mark(x_local,y_local,2)==1)
                j_1=out_fic_index(x_local+1,y_local);
                if(out_fic_vari_x(j_1,2)==0)
                    out_fic_right_x=out_fic_vari_x(j_1,1);
                    out_fic_right_y=out_fic_vari_y(j_1,1);
                elseif(out_fic_vari_x(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==2)
                            out_fic_right_x=out_fic_vari_x(j_1,k);
                            out_fic_right_y=out_fic_vari_y(j_1,k);
                            break;
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+d_d(x_global+1,y_global,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_right_y-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+d_d(x_global+1,y_global,2)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_right_x-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            elseif(in_cross_boun_mark(x_local,y_local,2)==2)
                j_1=out_fic_index(x_local,y_local+1);
                if(out_fic_vari_x(j_1,2)==0)
                    out_fic_upper_x=out_fic_vari_x(j_1,1);
                    out_fic_upper_y=out_fic_vari_y(j_1,1);
                elseif(out_fic_vari_x(j_1,2)~=0)
                    for k=1:1:2
                        if(out_dir_index(j_1,k)==1)
                            out_fic_upper_x=out_fic_vari_x(j_1,k);
                            out_fic_upper_y=out_fic_vari_y(j_1,k);
                            break;
                        end
                    end
                end
                equilibrium(x_global,y_global,1)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global-1,y_global,1)+out_fic_vari_x(i,1)-2*d_d(x_global,y_global,1))/(d_x^2)+...
                    (1-2*v_aluminum)*(d_d(x_global,y_global-1,1)+d_d(x_global,y_global+1,1)-2*d_d(x_global,y_global,1))/(d_y^2)+...
                    (out_fic_upper_y-d_d(x_global-1,y_global+1,2)-d_d(x_global+1,y_global-1,2)+d_d(x_global-1,y_global-1,2))/(4*d_x*d_y));
                equilibrium(x_global,y_global,2)=(miu_aluminum+lambda_aluminum)*(2*(1-v_aluminum)*...
                    (d_d(x_global,y_global-1,2)+d_d(x_global,y_global+1,2)-2*d_d(x_global,y_global,2))/(d_y^2)+...
                    (1-2*v_aluminum)*(d_d(x_global-1,y_global,2)+out_fic_vari_y(i,1)-2*d_d(x_global,y_global,2))/(d_x^2)+...
                    (out_fic_upper_x-d_d(x_global-1,y_global+1,1)-d_d(x_global+1,y_global-1,1)+d_d(x_global-1,y_global-1,1))/(4*d_x*d_y));
            end
        end
    end
    
    
    for i=1:1:1
        for j=1:1:max_y
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=0;
            d_d(i,j,1)=0;
            d_d(i,j,2)=0;
        end
    end
    for j=1:1:1
        for i=1:1:max_x
            equilibrium(i,j,1)=0;
            equilibrium(i,j,2)=0;
            d_d(i,j,1)=0;
            d_d(i,j,2)=0;
        end
    end
    
	equilibrium(max_x,1,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
		lambda_steel*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
    equilibrium(max_x,1,2)=lambda_steel*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
            (2*miu_steel+lambda_steel)*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
	
	equilibrium(max_x,max_y,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
		lambda_steel*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
	equilibrium(max_x,max_y,2)=lambda_steel*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
        (2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
		
	for j=2:1:max_y-1
		equilibrium(max_x,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,j,1)-4*d_d(max_x-1,j,1)+d_d(max_x-2,j,1))/(2*d_x)+...
			lambda_steel*(d_d(max_x,j+1,2)-d_d(max_x,j-1,2))/(2*d_y);
		equilibrium(max_x,j,2)=miu_steel*(d_d(max_x,j+1,1)-d_d(max_x,j-1,1))/(2*d_y)+...
			miu_steel*(3*d_d(max_x,j,2)-4*d_d(max_x-1,j,2)+d_d(max_x-2,j,2))/(2*d_x);
	end
	for i=2:1:max_x-1
		equilibrium(i,max_y,2)=lambda_steel*(d_d(i+1,max_y,1)-d_d(i-1,max_y,1))/(2*d_x)+...
			(2*miu_steel+lambda_steel)*(3*d_d(i,max_y,2)-4*d_d(i,max_y-1,2)+d_d(i,max_y-2,2))/(2*d_y);
		equilibrium(i,max_y,1)=miu_steel*(3*d_d(i,max_y,1)-4*d_d(i,max_y-1,1)+d_d(i,max_y-2,1))/(2*d_y)+...
			miu_steel*(d_d(i+1,max_y,2)-d_d(i-1,max_y,2))/(2*d_x);
	end
%{
	for i=2:1:high_x-2
		for j=1:1:1
			equilibrium(i,j,1)=miu_aluminum*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_aluminum*(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=lambda_aluminum*(d_d(i+1,1,1)-d_d(i-1,1,1))/(2*d_x)+...
                (2*miu_aluminum+lambda_aluminum)*(-3*d_d(i,1,2)+4*d_d(i,2,2)-d_d(i,3,2))/(2*d_y);
		end
	end
	for i=high_x-1:1:high_x-1
		for j=1:1:1
			equilibrium(i,j,1)=miu_aluminum*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_aluminum*(out_fic_vari_y(size_out_fic_for_in,1)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=lambda_aluminum*(out_fic_vari_x(size_out_fic_for_in,1)-d_d(i-1,1,1))/(2*d_x)+...
        (2*miu_aluminum+lambda_aluminum)*(-3*d_d(i,1,2)+4*d_d(i,2,2)-d_d(i,3,2))/(2*d_y);
		end
	end
	for i=high_x:1:high_x
		for j=1:1:1
			equilibrium(i,j,1)=miu_steel*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,j,2)-in_fic_vari_y(size_in_fic_for_out,1))/(2*d_x);
			equilibrium(i,j,2)=lambda_steel*(d_d(i+1,1,1)-in_fic_vari_x(size_in_fic_for_out,1))/(2*d_x)+...
        (2*miu_steel+lambda_steel)*(-3*d_d(i,1,2)+4*d_d(i,2,2)-d_d(i,3,2))/(2*d_y);
		end
	end
	for i=high_x+1:1:max_x-1
		for j=1:1:1
			equilibrium(i,j,1)=miu_steel*(-3*d_d(i,j,1)+4*d_d(i,j+1,1)-d_d(i,j+2,1))/(2*d_y)+...
                miu_steel*(d_d(i+1,j,2)-d_d(i-1,j,2))/(2*d_x);
			equilibrium(i,j,2)=lambda_steel*(d_d(i+1,1,1)-d_d(i-1,1,1))/(2*d_x)+...
        (2*miu_steel+lambda_steel)*(-3*d_d(i,1,2)+4*d_d(i,2,2)-d_d(i,3,2))/(2*d_y);
		end
    end
    
	equilibrium(max_x,1,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
		lambda_steel*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
    equilibrium(max_x,1,2)=lambda_steel*(3*d_d(max_x,1,1)-4*d_d(max_x-1,1,1)+d_d(max_x-2,1,1))/(2*d_x)+...
            (2*miu_steel+lambda_steel)*(-3*d_d(max_x,1,2)+4*d_d(max_x,2,2)-d_d(max_x,3,2))/(2*d_y);
	
	equilibrium(max_x,max_y,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
		lambda_steel*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
	equilibrium(max_x,max_y,2)=lambda_steel*(3*d_d(max_x,max_y,1)-4*d_d(max_x-1,max_y,1)+d_d(max_x-2,max_y,1))/(2*d_x)+...
        (2*miu_steel+lambda_steel)*(3*d_d(max_x,max_y,2)-4*d_d(max_x,max_y-1,2)+d_d(max_x,max_y-2,2))/(2*d_y);
		
	for j=2:1:max_y-1
		equilibrium(max_x,j,1)=(2*miu_steel+lambda_steel)*(3*d_d(max_x,j,1)-4*d_d(max_x-1,j,1)+d_d(max_x-2,j,1))/(2*d_x)+...
			lambda_steel*(d_d(max_x,j+1,2)-d_d(max_x,j-1,2))/(2*d_y);
		equilibrium(max_x,j,2)=miu_steel*(d_d(max_x,j+1,1)-d_d(max_x,j-1,1))/(2*d_y)+...
			miu_steel*(3*d_d(max_x,j,2)-4*d_d(max_x-1,j,2)+d_d(max_x-2,j,2))/(2*d_x);
	end
	for i=2:1:max_x-1
		equilibrium(i,max_y,2)=lambda_steel*(d_d(i+1,max_y,1)-d_d(i-1,max_y,1))/(2*d_x)+...
			(2*miu_steel+lambda_steel)*(3*d_d(i,max_y,2)-4*d_d(i,max_y-1,2)+d_d(i,max_y-2,2))/(2*d_y);
		equilibrium(i,max_y,1)=miu_steel*(3*d_d(i,max_y,1)-4*d_d(i,max_y-1,1)+d_d(i,max_y-2,1))/(2*d_y)+...
			miu_steel*(d_d(i+1,max_y,2)-d_d(i-1,max_y,2))/(2*d_x);
	end
	%}
    
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
    %{
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
	%}
    for i=1:1:max_x
        for j=1:1:max_y
            if((i==1)||(j==1))
                coefficient_displacement(i,j,1)=0;
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
	coefficient_b=zeros(max_x,max_y,2);
	for i=1:1:max_x
		for j=1:1:max_y
			if(i==max_x)
				coefficient_b(i,j,2)=0;%1000N force added on the boundary
				coefficient_b(i,j,1)=1000;
			else
				coefficient_b(i,j,1)=0;
				coefficient_b(i,j,2)=0;
			end
		end
	end

	total_displacement=sym(zeros(max_x,max_y,2));
	for i=1:1:max_x
		for j=1:1:max_y
			for k=1:1:2
				if(k==1)
					total_displacement(i,j,k)=Ux(i,j);
				else
					total_displacement(i,j,k)=Uy(i,j);
				end
			end
		end
	end

	%iterations
	accuracy=1e-15;
	displacement=zeros(max_x,max_y,2);
	updated_displacement=zeros(max_x,max_y,2);

	for i=1:1:max_x
		for j=1:1:max_y
			for k=1:1:2
				if(coefficient_displacement(i,j,k)==0)
					displacement(i,j,k)=0;
				end
			end
		end
	end

	%method to parse variables in the equilibrium before enter into the
	%iteration process

	variable_coefficient_cell=cell(max_x,max_y,2);
for i=1:1:max_x
    for j=1:1:max_y
        for k=1:1:2
            variables=symbolic_variable_cell{i,j,k};
            coeff_variables=zeros(length(variables),1);
            %get coefficient of each variables in each equilibrium
            for m=1:1:length(variables)
                t=coeffs(equilibrium(i,j,k),variables(m));
                if(size(t,2)==1)
                    coeff_variables(m)=t;
                else
                    coeff_variables(m)=t(2);
                end
            end
            variable_coefficient_cell{i,j,k}=coeff_variables;
        end
    end
end

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
                        updated_displacement(i,j,k)=displacement(i,j,k)-1.15*(residue-coefficient_b(i,j,k))/coefficient_displacement(i,j,k);
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
end