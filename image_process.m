function[]=image_process()
    colormap jet;
    t=colorbar;
    k=get(t,'Limits');
    set(t,'Ticks',linspace(k(1),k(2),10));
    xlabel('x position');
    ylabel('y position');
    zlabel('stress');
end