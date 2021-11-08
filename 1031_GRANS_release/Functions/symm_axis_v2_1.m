function ax_lim=symm_axis_v2_1(f,perc)
    B=prctile(abs(vert(f)),perc);
    B=10^ceil(log10(B));
    ax_lim=B*[-1 1]; 
end
