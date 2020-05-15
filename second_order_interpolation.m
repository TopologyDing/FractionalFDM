function[results]=second_order_interpolation(disp,posi,center_coord,x_or_y,d_x,d_y,options)
    syms a b c;
    %option:
    %option==1: compute displacement
    %option==2: compute derivative
    %option==3: compute both
    %x_or_y:
    %x_or_y==1: compute x direction
    %x_or_y==2: compute y direction
    if(x_or_y==1)% x direction
        eqn(1)=a*((center_coord-2)*d_x)^2+b*((center_coord-2)*d_x)+c==disp(1);
        eqn(2)=a*((center_coord-1)*d_x)^2+b*((center_coord-1)*d_x)+c==disp(2);
        eqn(3)=a*((center_coord)*d_x)^2+b*((center_coord)*d_x)+c==disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        if(options==1)
            results=a*posi^2+b*posi+c;
        elseif(options==2)
            results=2*a*posi+b;
        elseif(options==3)
            results(1)=a*posi^2+b*posi+c;
            results(2)=2*a*posi+b;
        end
    else % y direction 
        eqn(1)=a*((center_coord-2)*d_y)^2+b*((center_coord-2)*d_y)+c==disp(1);
        eqn(2)=a*((center_coord-1)*d_y)^2+b*((center_coord-1)*d_y)+c==disp(2);
        eqn(3)=a*((center_coord)*d_y)^2+b*((center_coord)*d_y)+c==disp(3);
        [a,b,c]=solve(eqn,[a,b,c]);
        if(options==1)
            results=a*posi^2+b*posi+c;
        elseif(options==2)
            results=2*a*posi+b;
        elseif(options==3)
            results(1)=a*posi^2+b*posi+c;
            results(2)=2*a*posi+b;
        end
    end
end