__author__ = 'Michael.Marty'


import mdtraj as md
import numpy as np
import os
import matplotlib.pyplot as plt
from bisect import bisect_left
import matplotlib.cm as cm
import matplotlib.gridspec as gridspec

def MakeCmap(tup,alpha):
    cdict = {'red':  ((0.0,tup[0], tup[0]),
                       (1.0, tup[0], tup[0])),
             'green': ((0.0,tup[1], tup[1]),
                       (1.0, tup[1], tup[1])),
             'blue':  ((0.0, tup[2], tup[2]),
                       (1.0, tup[2], tup[2]))}
    cdict['alpha'] = ((0.0, 0, 0),
                       (1.0, alpha, alpha))
    return cdict

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

def nearestunsorted(A,target):
    return np.argmin(np.abs(A-target))

def nearest(A, target):
    i = bisect_left(A, target)
    if i<=0:
        return 0
    elif i>=len(A)-1:
        return len(A)-1
    if np.abs(A[i] - target) > np.abs(A[i-1]-target):
        i-=1
    return i

def conv2d(A,B,norm=False):
    a=np.fft.fft2(A)
    b=np.fft.fft2(B)
    c=a*b
    out=np.fft.ifft2(c).real
    if norm:
        out=(np.clip(out,0.1,100000)-0.1)
    return out

def smarthist(a,binsize):
    ints=np.round(a/binsize)
    bincount=np.bincount(ints.astype(int))
    xvals=np.arange(0,len(bincount))*binsize
    return np.transpose([xvals,bincount,np.cumsum(bincount)])

def classify(data,hpos,ppos,npos,lpos,spos):
    types=[]
    for i,row in enumerate(data):
        columns=[3,10,12,2,5]
        cutoffs=[2.2,4,3.25,2.1,10]
        type=[0,0,0,0,0,0,0]
        for j,c in enumerate(columns):
            if row[c]<cutoffs[j]:
                if j<3:
                    type[j]=1
                else:
                    type[j+2]=1
        if type[1]==1 or type[2]==1:
            type[0]=0
            type[3]=1
            type[4]=1
        if type[0]==1:
            type[4]=1
        if type[5]==1:
            type[6]=0
        types.append(type)

        if type[0]==1:
            hpos.append([row[6],row[7],row[8]])
        if type[1]==1:
            ppos.append([row[6],row[7],row[8]])
        if type[2]==1:
            npos.append([row[6],row[7],row[8]])
        if type[5]==1:
            lpos.append([row[6],row[7],row[8]])
        if type[6]==1:
            spos.append([row[6],row[7],row[8]])

    types=np.array(types)
    out=np.sum(types,axis=0)
    return out

#Types:
#0: H-Bond
#1: LPos
#2: LNeg
#3: Ionic
#4: Ionic+Hbond
#5: Shell
#6: Second Shell


#0: Lipid Number
#1: Global Min Distance All Lipid - All Protein
#2: Min Distance Lipid H - Prot H
#3: Min Distance Lipid NoH - Prot H
#4: Min Distance Lipid H - Prot NoH
#5: Min Distance Lipid NoH - Prot NoH
#6 7 8: Lipid Average X Y Z
#9: Min Distance Lipid NTL - Prot Neg
#10: Min Distance Lipid Choline Methyl Carbons - Prot Neg
#11: Min Distance Lipid P - Prot Pos
#12: Min Distance Lipid Phosphate Oxygens - Prot Pos

terminals=[["(residue 1 and name N)","(residue 231 and name C)"],["(residue 2 and name N)","(residue 386 and name C)"]]

labels=["Lipid Index","Any Atoms","HH","NoH-H","H-NoH","NoH-NoH","X","Y","Z","LN+Pneg","LMethylC-Pneg","LP-Ppos","LPO-Ppos"]

proteins=["AqpZ","AmtB"]
#Note: AqpZ is oriented with the top as cytoplasm, flip it so the bottom is cytoplasm
inversion=[-1,1]


starts=[200,200]
ends=[1020,1019]

ncol=13

All=[]
AllSep=[]
for c in range(0,ncol):
    All.append([])

gs=gridspec.GridSpec(2,2)
gs.update(wspace=0.05,hspace=0.15)
fig=plt.figure(figsize=(4.13,4.5))
size=8

for i,protein in enumerate(proteins):
    if i==i:
        classes=[]
        hpos=[]
        ppos=[]
        npos=[]
        lpos=[]
        spos=[]
        All2=[]
        for c in range(0,ncol):
            All2.append([])

        frames=np.arange(starts[i],ends[i])
        dir="C:\\VMD\\"+protein+"\\AnalysisAlignedLong"
        os.chdir(dir)
        for j,f in enumerate(frames):

            file="HSalt"+str(f)+".txt"
            path=os.path.join(dir,file)
            data=np.loadtxt(path)
            data[:,8]=data[:,8]*inversion[i]

            for c in range(0,ncol):
                dat=data[data[:,c]<10,c]
                All[c]=np.hstack((All[c],dat))
                All2[c]=np.hstack((All2[c],dat))


            classes.append(classify(data,hpos,ppos,npos,lpos,spos))
        times=frames*0.1
        AllSep.append(All2)

        classes=np.array(classes)
        cavg=np.average(classes,axis=0)
        cstd=np.std(classes,axis=0)
        #print cavg
        #print cstd


        dir="C:\\VMD\\"+protein
        os.chdir(dir)
        print("Loading...")
        t=md.load_dcd("proteinaligned.dcd",top="protein.psf")


        sel=t.topology.select("all")
        back=t.atom_slice(sel,inplace=False)

        sel=t.topology.select(terminals[i][1]+"or ((resname GLU or resname ASP) and (name CD or name CG))")
        neg=t.atom_slice(sel)

        sel=t.topology.select(terminals[i][0]+"or ((resname LYS or resname ARG) and (name NZ or name CZ))")
        pos=t.atom_slice(sel)
        '''
        ncoords=neg.xyz[200:]*10.
        navg=np.average(ncoords,axis=0)
        ntop=navg[navg[:,2]>0]
        nbottom=navg[navg[:,2]<0]
        nparts=[ntop,nbottom]


        pcoords=pos.xyz[200:]*10.
        pavg=np.average(pcoords,axis=0)
        ptop=pavg[pavg[:,2]>0]
        pbottom=pavg[pavg[:,2]<0]
        pparts=[ptop,pbottom]
        print navg.shape,pavg.shape
        '''

        coords=back.xyz[starts[i]:]*10.
        coords[:,2]=coords[:,2]*inversion[i]
        avg=np.average(coords,axis=0)
        atop=avg[avg[:,2]>0]
        abottom=avg[avg[:,2]<0]

        minmax=[-60,60]
        xaxis=np.arange(minmax[0],minmax[1],step=0.1)
        yaxis=np.arange(minmax[0],minmax[1],step=0.1)
        xgrid,ygrid=np.meshgrid(xaxis,yaxis,indexing="ij")
        rgrid,pgrid=cart2pol(xgrid,ygrid)

        kgrid=np.zeros_like(xgrid)
        for k in [0,1]:
            for j in [0,1]:
                center=[minmax[k],minmax[j]]
                bool1=np.sqrt((xgrid-center[0])**2+(ygrid-center[1])**2)<1.5
                kgrid[bool1]=1

        step=2*np.pi/360.
        paxis=np.arange(-np.pi,np.pi,step=step)

        parts=[atop,abottom]
        traces=[]
        pgrids=[]
        ngrids=[]
        for k,avgpart in enumerate(parts):
            igrid=np.zeros_like(xgrid)
            posgrid=np.zeros_like(xgrid)
            neggrid=np.zeros_like(xgrid)
            for coor in avgpart:
                x=nearest(xaxis,coor[0])
                y=nearest(yaxis,coor[1])
                igrid[x,y]=1
            '''
            for coor in pparts[i]:
                x=nearest(xaxis,coor[0])
                y=nearest(yaxis,coor[1])
                posgrid[x,y]=1

            for coor in nparts[i]:
                x=nearest(xaxis,coor[0])
                y=nearest(yaxis,coor[1])
                neggrid[x,y]=1
            neggrid=conv2d(neggrid,kgrid,norm=True)
            posgrid=conv2d(posgrid,kgrid,norm=True)

            '''

            print("Convolving...")
            igrid=conv2d(igrid,kgrid)

            bool1=igrid>0.5
            pgrid2=pgrid[bool1]
            rgrid2=rgrid[bool1]

            raxis=[]
            for p in paxis:
                bool2=np.abs(pgrid2-p)<step/2.
                maxr=np.amax(rgrid2[bool2])
                raxis.append(maxr)
            xtrace,ytrace=pol2cart(raxis,paxis)
            traces.append([xtrace,ytrace])

            for j in range(starts[i],ends[i]-1):
                ncoords=neg.xyz[j]*10.
                ncoords[:,2]=ncoords[:,2]*inversion[i]
                ntop=ncoords[ncoords[:,2]>0]
                nbottom=ncoords[ncoords[:,2]<0]
                nparts=[ntop,nbottom]

                pcoords=pos.xyz[j]*10.
                pcoords[:,2]=pcoords[:,2]*inversion[i]
                ptop=pcoords[pcoords[:,2]>0]
                pbottom=pcoords[pcoords[:,2]<0]
                pparts=[ptop,pbottom]

                for coor in pparts[k]:
                    x=nearest(xaxis,coor[0])
                    y=nearest(yaxis,coor[1])
                    posgrid[x,y]+=1

                for coor in nparts[k]:
                    x=nearest(xaxis,coor[0])
                    y=nearest(yaxis,coor[1])
                    neggrid[x,y]+=1

            print("Convolving...")

            neggrid=conv2d(neggrid,kgrid,norm=True)
            posgrid=conv2d(posgrid,kgrid,norm=True)
            power=0.5
            neggrid=(neggrid/np.amax(neggrid))**power
            posgrid=(posgrid/np.amax(posgrid))**power

            newalpha=0.5
            cm.register_cmap(name="reds",data=MakeCmap((1,0,0),newalpha))
            cm.register_cmap(name="blues",data=MakeCmap((0,0,1),newalpha))
            pgrids.append(posgrid)
            ngrids.append(neggrid)

        hpos=np.array(hpos)
        ppos=np.array(ppos)
        npos=np.array(npos)
        lpos=np.array(lpos)
        spos=np.array(spos)
        alpha=0.25
        alpha2=1

        axtop=fig.add_subplot(gs[i*2],aspect="equal")
        plt.setp(axtop.get_xticklabels(),fontsize=size)
        plt.setp(axtop.get_yticklabels(),fontsize=size)

        lw=2
        linecol="k"
        ls="-"
        xtrace,ytrace=traces[0]
        zo=0
        plt.plot(xtrace,ytrace,color=linecol,linewidth=lw,linestyle=ls,zorder=zo)

        plt.contourf(xgrid,ygrid,pgrids[0],100,cmap="blues")
        plt.contourf(xgrid,ygrid,ngrids[0],100,cmap="reds")
        boo1=hpos[:,2]>0
        boo2=ppos[:,2]>0
        boo3=npos[:,2]>0
        boo4=lpos[:,2]>0
        boo5=spos[:,2]>0
        #plt.scatter(spos[boo5,0],spos[boo5,1],c="w",alpha=alpha/4.)

        #calculate the scatter size based on the area of a lipid
        diffs=(axtop.transData.transform((np.sqrt(59/np.pi),0))-axtop.transData.transform((0,0)))
        s=diffs[0]**2

        plt.scatter(lpos[boo4,0],lpos[boo4,1],c="w",s=s,alpha=alpha/4.)

        plt.scatter(hpos[boo1,0],hpos[boo1,1],c="w",s=s,alpha=alpha/4.)
        plt.scatter(ppos[boo2,0],ppos[boo2,1],c="r",s=s,alpha=alpha)
        plt.scatter(npos[boo3,0],npos[boo3,1],c="b",s=s,alpha=alpha)

        periplasm=np.array([len(lpos[boo4])+len(hpos[boo1]),len(ppos[boo2]),len(npos[boo3])])/float(len(frames))

        #plt.scatter(atop[:,0],atop[:,1],c="k",alpha=alpha2)
        #plt.scatter(ptop[:,0],ptop[:,1],c="b",alpha=alpha2)
        #plt.scatter(ntop[:,0],ntop[:,1],c="r",alpha=alpha2)

        plt.title(protein + " Periplasm",fontsize=size)
        if i>0:
            plt.xlabel("X Position (${\AA}$)",fontsize=size)
        else:
            axtop.set_xticklabels([])

        plt.ylabel("Y Position (${\AA}$)",fontsize=size)
        plt.gca().set_ylim(-60,60)
        plt.gca().set_xlim(-60,60)


        axbottom=fig.add_subplot(gs[i*2+1],aspect="equal")
        plt.setp(axbottom.get_xticklabels(),fontsize=size)
        plt.setp(axbottom.get_yticklabels(),fontsize=size)
        plt.gca().set_ylim(-60,60)
        plt.gca().set_xlim(-60,60)

        if i>0:
            plt.xlabel("X Position (${\AA}$)",fontsize=size)
            axbottom.set_xticklabels(np.arange(-40,80,20))
        else:
            axbottom.set_xticklabels([])

        axbottom.set_yticklabels([])



        plt.contourf(xgrid,ygrid,pgrids[1],100,cmap="blues")
        plt.contourf(xgrid,ygrid,ngrids[1],100,cmap="reds")

        boo1=hpos[:,2]<0
        boo2=ppos[:,2]<0
        boo3=npos[:,2]<0
        boo4=lpos[:,2]<0
        boo5=spos[:,2]<0
        #plt.scatter(spos[boo5,0],spos[boo5,1],c="w",alpha=alpha/4.)
        plt.scatter(lpos[boo4,0],lpos[boo4,1],c="w",s=s,alpha=alpha/4.)

        plt.scatter(hpos[boo1,0],hpos[boo1,1],c="w",s=s,alpha=alpha/4.)
        plt.scatter(ppos[boo2,0],ppos[boo2,1],c="r",s=s,alpha=alpha)
        plt.scatter(npos[boo3,0],npos[boo3,1],c="b",s=s,alpha=alpha)

        cytoplasm=np.array([len(lpos[boo4])+len(hpos[boo1]),len(ppos[boo2]),len(npos[boo3])])/float(len(frames))

        #plt.scatter(abottom[:,0],abottom[:,1],c="k",alpha=alpha2)
        #plt.scatter(pbottom[:,0],pbottom[:,1],c="b",alpha=alpha2)
        #plt.scatter(nbottom[:,0],nbottom[:,1],c="r",alpha=alpha2)

        xtrace,ytrace=traces[1]
        plt.plot(xtrace,ytrace,color=linecol,linewidth=lw,linestyle=ls,zorder=zo)

        plt.title(protein + " Cytoplasm",fontsize=size)


        print(periplasm)
        print(cytoplasm)
        print(periplasm/(periplasm+cytoplasm))
        print(periplasm+cytoplasm)

plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_4.png")
#plt.savefig("C:\\Users\\michael.marty\\Dropbox\\Research\\MSNanodisc\\Figures\\Figure_4.eps")
plt.show()




