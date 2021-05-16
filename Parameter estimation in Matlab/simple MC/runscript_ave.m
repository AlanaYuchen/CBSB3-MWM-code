% start run estpar
bestres = zeros(4,8,4,20);
load allvars
load alldata

alltres = cell(1,3);
t=zeros(1,4);
mno = 0;
    for   strain = 1:2
        for temp = 1:2
                tic
                toc
                t1=toc;
            mno = mno + 1;
            disp(['Mouse S',num2str(strain),'T',num2str(temp)]); % convert numeric to string
            
%           [Wmult, sigma_pc, sigma_ac, Vdecay, ac_const, beta, etdecay, lrate, discf, noise, gof, PMs] = estpar (mousedata);
            [best_Wmult, best_sigma_pc, best_sigma_ac, best_Vdecay, ...
            best_ac_const, best_beta, best_etdecay, best_alpha, ...
            best_gamma, best_noise, best_wpun, best_gof, best_PMs, goferr_set, params, PMs]...
            = estpar_iter3 (squeeze(mean(squeeze(alldata(:,strain,temp,:,:,:)),2)),...
            allvars,strain,temp);
        
            bestres(mno,:,:,1) = reshape(best_Wmult,1,8,4,1);
            bestres(mno,:,:,2) = reshape(best_sigma_pc,1,8,4,1);
            bestres(mno,:,:,3) = reshape(best_sigma_ac,1,8,4,1);
            bestres(mno,:,:,4) = reshape(best_Vdecay,1,8,4,1);
            bestres(mno,:,:,5) = reshape(best_ac_const,1,8,4,1);
            bestres(mno,:,:,6) = reshape(best_beta,1,8,4,1);
            bestres(mno,:,:,7) = reshape(best_etdecay,1,8,4,1);
            bestres(mno,:,:,8) = reshape(best_alpha,1,8,4,1);
            bestres(mno,:,:,9) = reshape(best_gamma,1,8,4,1);
            bestres(mno,:,:,10) = reshape(best_noise,1,8,4,1);
            bestres(mno,:,:,11) = reshape(best_wpun,1,8,4,1);
            bestres(mno,:,:,12) = reshape(best_gof,1,8,4,1);
            bestres(mno,:,:,13:20) = reshape(best_PMs,1,8,4,8);
            
            alltres{mno, 1} = goferr_set;
            alltres{mno, 2} = params;
            alltres{mno, 3} = PMs; 
            toc
            t2=toc;
            t = t2-t1
            timeset(1,mno)= t;
            file = ['Group',num2str(mno),'.mat'];
            save(file)
        end
     end 
 
    
    save alltressimpleMC alltres
    save bestressimpleMC bestres

            