        d_dens = dens_00_data - dens_data

        #plot results
        fig, ax = plt.subplot_mosaic([['top','top'],['left', 'right']],
                                    gridspec_kw={'height_ratios': [1,3]})

        ax['top'].plot(tgrid*au2fs,Efield*au2Vnm)
        ax['top'].set_xlim([tmin*au2fs,tmax*au2fs])
        ax['top'].set_xlabel("Time [fs]",labelpad=10)
        ax['top'].set_ylabel("E [V/nm]")

        time = tgrid[it]*au2fs
        time_str = "Time: %s fs" % ('{:6.2f}'.format(time))

        plt.figtext(0.75,0.93,time_str)

        ax['top'].axvline(tgrid[it]*au2fs,c='r',lw=3)

        ax['left'].contourf(kx_grid/au2nm,ky_grid/au2nm,rho_data,levels=np.linspace(0.,1.,30),cmap=cm.jet)
        ax['left'].set_box_aspect(1)
        ax['left'].set_xlabel(r"$k_x [\mathrm{nm}^{-1}]$",labelpad=10)
        ax['left'].set_ylabel(r"$k_y [\mathrm{nm}^{-1}]$",labelpad=5)

        ax['right'].contourf(xgrid*au2A,ygrid*au2A,d_dens,levels=np.linspace(-2e-12,2e-12,30),cmap=cm.bwr)
        ax['right'].set_box_aspect(1)
        ax['right'].set_xlabel(r"$x [\AA]$",labelpad=10)
        ax['right'].set_ylabel(r"$y [\AA]$",labelpad=-3)   
        ax['right'].set_xlim([xmin*au2A,xmax*au2A])
        ax['right'].set_ylim([ymin*au2A,ymax*au2A])

        for posA in np.concatenate((A1,A2)):
            #ax['right'].plot(posA[0],posA[1],'ro')
            ax['right'].scatter(posA[0],posA[1],marker="o",s=70,alpha=0.4,c='black')

        #for posA in A2:
        #    ax['right'].plot(posA[0],posA[1],'ro')        

        plt.subplots_adjust(left=0.1, bottom=0.04, right=0.9, top=0.97, wspace=0.2, hspace=0.1)

        plt.savefig(out_t_fname,dpi=150)
        plt.close()