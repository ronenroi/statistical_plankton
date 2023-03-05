close all
w=1;
xbins3 = -90:w:90;

for i = 1:length(n_theta_vec)
    figure
    for j = 1:length(max_scales)
        for k = 1:length(max_shifts)
            n_theta = n_theta_vec(i);

            [thetadiff ]=comparecl_model( clstack(:,:,k,j,i), ...
                ref_clmatrix, n_theta );
            thetadiff(logical(eye(size(thetadiff)))) = []; % Or A = A(~eye(size(A)))
            thetadiff=thetadiff-360/n_theta_vec(i)+1;
            N_lines = length(thetadiff);
            
            f = @(par,x)(N_lines * par(1))* w  * 1/(par(2)*sqrt(2*pi)) * (exp(-x.^2 ./(2*par(2)^2))+exp(-(180-x).^2./(2*par(2)^2))) ...
                + (N_lines * (1- par(1)))* w  / 180;
            B0 = [0.8, 1.4];
            opts = statset('nlinfit');
            opts.RobustWgtFun = 'cauchy';
            h=hist(thetadiff,xbins3);
            B = nlinfit(xbins3(1:180), h(1:180), f, B0);
            
            subplot(length(max_scales),length(max_shifts),k+(j-1)*length(max_shifts));hist(thetadiff,xbins3)
            hold on
            plot(xbins3, f(B,xbins3), '-r', 'LineWidth',1.5)
            hold off
            legend(sprintf('P %f sigma %f',B(1),B(2)));
            P(k,j,i)=B(1);
            sigma(k,j,i) = B(2);
            %     e= h(2:181) - f(B,-89:90);
            %     ss_tot = sum((h - mean(h)).^2);
            %     R_vec(i,j)= 1 - sum(e.^2)/ss_tot;
            
        end
    end
end