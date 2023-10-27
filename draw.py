import matplotlib.pyplot as plt
import numpy as np



def draw_cd(pval, avg_rank, barwidth=0.75,linew=2,deltasep=1.0,deltasame=0.2,lineoffs=0.3,cuffmult=0.25,labels_inside=True,label_size=20,link_cuffs=0.2,textmult=0.4,lowerx=0.8,gridx=0.2,additional_text=None):
    """expects a list of comparisons (algo1,algo2,pval,significant), and a dictionary for the average ranks
    barwidth: width of the bars
    linew: width of the lines
    deltasep: when not same group, how big is the gap
    deltasame: when same group, how big is the gap
    lineoffs: how far from the bar is the line to show connectivity
    cuffmult: how big is the cuff
    labels_inside: whether to put the labels inside the bars or just use as ticks
    label_size: size of the labels if drawn inside
    link_cuffs: wether to draw a small line between cuffs and bars. Either alpha value or False to not draw them
    textmult: how far away is the text, only matters if labels_inside
    lowerx: start of the x axis
    gridx: draw vertical lines at the integer ranks, boolean or alpha value
    additional_text: add these texts to each name. Either None, String or dictionary

    """
    algos=list(avg_rank.keys())
    ranks=[avg_rank[algo] for algo in algos]
    flipper=int(np.ceil(np.max(ranks)))+1
    ranks=[flipper-rank for rank in ranks]   
    #sort both
    ranks, algos = zip(*sorted(zip(ranks, algos),reverse=True))
    #horizontal bar plot

    if additional_text is None:
        additional_text={algo:"" for algo in algos}
    if type(additional_text)==str:
        additional_text={algo:additional_text for algo in algos}

    alg_to_ypos = {alg: ypos for ypos, alg in enumerate(algos)}
    #alg_to_upper = {alg: ypos + barwidth / 2 for alg, ypos in alg_to_ypos.items()}
    #alg_to_lower = {alg: ypos - barwidth / 2 for alg, ypos in alg_to_ypos.items()}

    alg_to_x= {alg: rank for alg, rank in zip(algos, ranks)}
    dex_to_algo={dex:alg for alg,dex in alg_to_ypos.items()}
    dex_to_x={dex:alg_to_x[dex_to_algo[dex]] for dex in dex_to_algo}





    #filter pval for pointless connections
    #need to reorder first
    todraw=[(alg1,alg2) for alg1,alg2,_,sign in pval if sign==False]
    todraw=[(alg1,alg2) if alg_to_ypos[alg1]>alg_to_ypos[alg2] else (alg2,alg1) for alg1,alg2 in todraw]

    basepos=[(alg_to_ypos[alg1],alg_to_ypos[alg2],alg1,alg2) for alg1,alg2 in todraw]

    def is_contained(a,b,arr):
        for mu,nu,_,_ in arr:
            if mu==a and nu==b:continue
            if mu>=a and nu<=b:return True
        return False

    seperators=[]
    for i in range(1,len(algos)):
        value=i-0.5
        if is_contained(value,value,basepos):continue
        #plt.axhline(y=value, color='black', linewidth=linew,linestyle='--')
        seperators.append(value)
    seperators=np.array(seperators)
    #print(seperators)
    #exit()

    alg_to_group={}
    for alg,ypos in alg_to_ypos.items():
        if len(seperators)==0:
            alg_to_group[alg]=0
        else:
            alg_to_group[alg]=np.sum(seperators<ypos)
    

    alg_to_upper = {}
    alg_to_lower = {}
    alg_to_center = {}
    loc=-barwidth/2
    deltasep+=barwidth
    deltasame+=barwidth
    last=None
    for alg,ypos in alg_to_ypos.items():
        if last is not None:
            if alg_to_group[alg]==alg_to_group[last]:
                loc+=deltasame
            else:
                loc+=deltasep

        alg_to_lower[alg]=loc
        alg_to_upper[alg]=loc+barwidth
        alg_to_center[alg]=loc+barwidth/2


        last=alg




    todraw=[(alg_to_center[alg1],alg_to_center[alg2],alg1,alg2) for alg1,alg2 in todraw]

    todraw=[(mu,nu,alg1,alg2) for mu,nu,alg1,alg2 in todraw if not mu==nu]
    todraw=[(mu,nu,alg1,alg2) if mu>nu else (nu,mu,alg2,alg1) for mu,nu,alg1,alg2 in todraw]
    todraw=[(mu,nu,alg1,alg2) for mu,nu,alg1,alg2 in todraw if not is_contained(mu,nu,todraw)]

    print(todraw)

    ypos=np.array([alg_to_lower[algo] for algo in algos])
    heis=np.array([alg_to_upper[algo]-alg_to_lower[algo] for algo in algos])

    xmax=max(ranks)
    xt=[i for i in range(1,int(np.ceil(xmax))+1)]# if i<=xmax]
    plt.xticks(xt,[flipper-zw for zw in xt])
    if gridx:
        for xx in xt:
            plt.axvline(x=xx, color='black', linewidth=1,linestyle='--',alpha=gridx,zorder=1)


    plt.barh(ypos, ranks, color='white', height=heis,align="edge")
    plt.barh(ypos, ranks, color='black', fill=False, linewidth=linew, edgecolor='black', height=heis,align="edge")

    maxx=np.min(ranks)
    count=len(todraw)
    delta=maxx/(count+1)


    def find_affected_dex(ymin,ymax):
        ymin=ymin-barwidth
        ymax=ymax+barwidth
        return [i for i,pos in enumerate(ypos) if pos>=ymin and pos<=ymax]

    #sort todraw so that: smaller lines first, then smaller x=inv smaller y first
    #first sort by average y
    todraw=sorted(todraw,key=lambda x: (x[0]+x[1])/2, reverse=True)
    #then sort by length
    todraw=sorted(todraw,key=lambda x: abs(x[0]-x[1]))


    start=delta
    curr=start
    baroffs=barwidth/2
    already_used={}

    for p1,p2,alg1,alg2 in todraw:
        if p1<p2:#so p1 is always >= p2
            p1,p2=p2,p1
            alg1,alg2=alg2,alg1
        affected=find_affected_dex(p2,p1)
        position=[dex_to_x[dex] for dex in affected]
        position=np.max(position)
        for dex in affected:
            if not dex in already_used:continue
            if np.max(already_used[dex])>position:
                position=np.max(already_used[dex])
        position+=lineoffs
        for dex in affected:
            if not dex in already_used:
                already_used[dex]=[]
            already_used[dex].append(position)
        plt.plot([position,position],[p1,p2], color='black', linewidth=linew)
        #little cuffs at the end
        plt.plot([position-cuffmult*lineoffs,position],[p1,p1], color='black', linewidth=linew)
        plt.plot([position-cuffmult*lineoffs,position],[p2,p2], color='black', linewidth=linew)
        #connect the cuffs to the xvalue
        if link_cuffs:
            plt.plot([alg_to_x[alg1],position-cuffmult*lineoffs],[p1,p1],color="black",linestyle="--",alpha=link_cuffs)
            plt.plot([alg_to_x[alg2],position-cuffmult*lineoffs],[p2,p2],color="black",linestyle="--",alpha=link_cuffs)
        curr+=delta

    #if two bars are not connected, draw an seperator line in between

    plt.xlabel("Average Rank")


    #mirror x axis both at the bottom and the top
    plt.tick_params(labeltop=True)

    if labels_inside:
        plt.tick_params(labelleft=False)
        for alg in algos:
            plt.text(alg_to_x[alg]-textmult*lineoffs,alg_to_center[alg],alg+additional_text[alg],ha="right",va="center",fontsize=label_size)

    else:
        yticks=[(alg_to_upper[algo]+alg_to_lower[algo])/2 for algo in algos]
        ylabels=[algo for algo in algos]
        plt.yticks(yticks, ylabels)


    plt.xlim(left=lowerx)

