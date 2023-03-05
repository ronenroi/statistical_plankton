%%
close all
w=1;
xbins3 = -90:w:90;
[ref_clmatrix,clcorr,cltheta]=clmatrix_cheat(data.inv_rot_matrices,360);
[clmatrix,clcorr,cltheta]=clmatrix_cheat(aligned_rots,360);

  [ thetadiff ]=comparecl_model( clmatrix, ...
        ref_clmatrix, 360);
    thetadiff(logical(eye(size(thetadiff)))) = []; % Or A = A(~eye(size(A)))
    thetadiff=thetadiff;
    %%%%%%%%%%%%%%
    N_lines = length(thetadiff);

    f = @(par,x)(N_lines * par(1))* w  * 1/(par(2)*sqrt(2*pi)) * (exp(-x.^2 ./(2*par(2)^2))+exp(-(180-x).^2./(2*par(2)^2))) ...
        + (N_lines * (1- par(1)))* w  / 180;
B0 = [0.8, 1.4];
 opts = statset('nlinfit');
 opts.RobustWgtFun = 'cauchy';
h=hist(thetadiff,xbins3);
B = nlinfit(xbins3(1:180), h(1:180), f, B0);

    hist(thetadiff,xbins3)
    hold on
    plot(xbins3, f(B,xbins3), '-r', 'LineWidth',1.5)
    hold off
    legend(sprintf('P %f sigma %f',B(1),B(2)));
    P=B(1);
    sigma= B(2);
    e= h(2:181) - f(B,-89:90);
    ss_tot = sum((h - mean(h)).^2);
    R_vec= 1 - sum(e.^2)/ss_tot;
        