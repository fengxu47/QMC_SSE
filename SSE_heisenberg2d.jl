#=
    SSE for 2d Heisenberg model(julia version)
    naive julia version of Anders Sandvik's ssebasic.f90
    refrence: http://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.html
              http://physics.bu.edu/~sandvik/programs/ssebasic/ssebasic.f90 
              http://physics.bu.edu/~sandvik/trieste12/tut1.pdf
              http://physics.bu.edu/~sandvik/trieste12/ssebasic.f90
               http://physics.bu.edu/~sandvik/trieste12/sseres.f90
=#

"""
bsites: label of two sites(spins) which is  connected by bond b

"""

function makelattice(nn::Int64,lx::Int64,ly::Int64,bsites::Matrix{Int64})
       
    for y1=0:ly-1
        for x1=0:lx-1
            s::Int64=1+x1+y1*lx
            x2=mod(x1+1,lx)
            y2=y1         #(x2,y2) right of (x1,y1)
            bsites[1,s]=s
            bsites[2,s]=1+x2+y2*lx  #horizontal bond
            x2=x1
            y2=mod(y1+1,ly)
            bsites[1,s+nn]=s
            bsites[2,s+nn]=1+x2+y2*lx #vertical bond
        end
    end
    return bsites
end


"""
diagonal update
1)add diagonal operator 
2)remove diagonal operator
3)off-diagonal operator=>flip spins
return nh

"""

function diagonalupdate(nh::Int64,mm::Int64,opstring::Vector{Int64},spin::Vector{Int64},bsites::Matrix{Int64},nb::Int64,beta::Int64)
    
    aprob :: Float64=0.50*beta*nb  #add diagonal operator prob
    dprob :: Float64=1.0/(0.5*beta*nb) #remove diagonal operator prob

    for i=1:mm
        op::Int64=opstring[i]
        # add diagonal operator
        if op==0 
            b = Int(trunc(nb*rand()))+1
            if spin[bsites[1,b]] != spin[bsites[2,b]]
                if (aprob>=(mm-nh) || aprob>=rand()*(mm-nh))
                    opstring[i]=2*b
                    nh=nh+1
                end
            end
        # remove diagonal operator
        elseif mod(op,2)==0
            p::Float64=dprob*(mm-nh+1)
            if (p>=1.0 || p>=rand())
                opstring[i]=0
                nh=nh-1
            end
        else
            b::Int64=fld(op,2)
            spin[bsites[1,b]]=-spin[bsites[1,b]]
            spin[bsites[2,b]]=-spin[bsites[2,b]]
        end
    end
    return nh
end


"""
construct linkvertex and loop update

"""
function loopupdate(mm::Int64,nn::Int64,opstring::Vector{Int64},bsites::Matrix{Int64},vertexlist::Vector{Int64},spin::Vector{Int64})
    

    firstspinop::Vector{Int64}=fill(-1,nn)
    lastspinop::Vector{Int64}=fill(-1,nn)

    #modify code of sandvik ssebasic.f90
    #main change: vertex legs's labels are 1, 2, 3, 4 rather than 0, 1, 2, 3

    #=
    construct linked vertex list
    =#
    for v0 in range(start=1,step=4,stop=4*mm-3)
        op=opstring[cld(v0,4)]
        if op != 0 # exist operator
            b=fld(op,2)
            s1=bsites[1,b]
            s2=bsites[2,b]
            v1=lastspinop[s1]
            v2=lastspinop[s2]
            if v1 != -1 #left(down) side of bond b is available for constructing a loop
                vertexlist[v1]=v0
                vertexlist[v0]=v1
            else
                firstspinop[s1]=v0
            end
            if v2 != -1 #right(up) side of bond b is avaaiable for constructing a loop
                vertexlist[v2]=v0+1
                vertexlist[v0+1]=v2
            else
                firstspinop[s2]=v0+1
            end
            lastspinop[s1]=v0+2
            lastspinop[s2]=v0+3
        else
            vertexlist[v0:v0+3]=[0, 0, 0, 0] #no operator assign vertexlist 0
        end
    end
    #creating the last links across the "time" boundary
    for s1=1:nn
        v1=firstspinop[s1]
        if v1 != -1
            v2=lastspinop[s1]
            vertexlist[v2]=v1
            vertexlist[v1]=v2
        end
    end


    #= 
    sweep of loop updates
    =#

    for v0 in range(start=1,step=2,stop=4*mm-1)
        if vertexlist[v0]<1  #remember we assian no operator's vertexlist zero
            continue
        end
        v1=v0
        if rand()<0.5
            while true
                opstring[cld(v1,4)]=xor(opstring[cld(v1,4)],1)    #diagonal operator <=> off-diagonal operator
                vertexlist[v1]=-1
                v2=xor(v1-1,1)+1 #from vertex v1 jump to its neighbor v2   vertex legs (1,2,3,4)
                v1=vertexlist[v2]
                vertexlist[v2]=-1 #-1=> the loop will be flipped
                if v1==v0
                    break
                end
            end
        else
            while true
                vertexlist[v1]=0
                v2=xor(v1-1,1)+1 #from vertex v1 jump to its neighbor v2   vertex legs (1,2,3,4)
                v1=vertexlist[v2]
                vertexlist[v2]=0
                if v1==v0
                    break
                end
            end
        end
    end


     for i=1:nn
        if firstspinop[i] != -1
            if vertexlist[firstspinop[i]]==-1  
                spin[i]=-spin[i] #-1: flip the spins in the loop
            end
        else
            if rand()<0.5 #a single spin no operator act on, single flip update
                spin[i]=-spin[i]
            end
        end
    end
    
end  #end function loopupdate


function measureobservables(mm::Int64,nh::Int64,nn::Int64,lx::Int64,spin::Vector{Int64},bsites::Matrix{Int64},opstring::Vector{Int64},enrg1::Float64,enrg2::Float64,amag2::Float64,ususc::Float64)


    am::Int64=0

    for i=1:nn
        am=am+spin[i]*((-1)^(mod(i-1,lx)+fld(i-1,lx))) #initial staggered magnetization
    end
    am=fld(am,2) 
    
    am2::Float64=0.0
    for i=1:mm
        op=opstring[i]
        if mod(op,2)==1 #off-diagonal operator
            b=fld(op,2)
            s1=bsites[1,b]
            s2=bsites[2,b]
            spin[s1]=-spin[s1]
            spin[s2]=-spin[s2]
            am=am+2*spin[s1]*((-1)^(mod(s1-1,lx)+fld(s1-1,lx))) 
            #off-diagonal operator flip a couple of spins,recalcuate staggered magnetization  
        end
        am2=am2+am^2

    end
    
    am2=am2/mm
    enrg1=enrg1+nh
    enrg2=enrg2+nh^2
    amag2=amag2+am2
    ususc=ususc+(sum(spin)/2)^2
    return enrg1,enrg2,amag2,ususc
end

function writeresults(nn::Int64,nb::Int64,msteps::Int64,beta::Int64,enrg1::Float64,enrg2::Float64,amag2::Float64,ususc::Float64)
    
    enrg1=enrg1/msteps
    enrg2=enrg2/msteps
    amag2=amag2/msteps
    ususc=ususc/msteps

    enrg2=(enrg2-enrg1*(enrg1+1.0))/nn
    enrg1=-(enrg1/(beta*nn)-0.25*nb/nn)
    amag2=3.0*amag2/(nn^2)
    ususc=beta*ususc/nn
    
    
    f=open("res.dat","a")
        println(f,enrg1,"   ",enrg2,"   ",amag2,"   ",ususc)
    close(f)
    # after write results, we should assign enrg1,enrg2,amag2,ususc zeros

end

"""
adjustcutoff: adjust mm

return mm
"""
function adjustcutoff(nh::Int64,mm::Int64,step::Int64,opstring::Vector{Int64},vertexlist::Vector{Int64})


    mmnew::Int64=nh+fld(nh,3)
    if mmnew>mm
        
        resize!(opstring,mmnew)
        opstring[mm+1:mmnew]=zeros(Int64,mmnew-mm)
        mm=mmnew
        resize!(vertexlist,4*mmnew)   

        f=open("cut.text","a")
            println(f,"step: ",step," Cut-off L: ",mmnew)
        close(f)
    end
    return mm
end


function mainprogram(lx::Int64,ly::Int64,beta::Int64,nbins::Int64,msteps::Int64,isteps::Int64)
    
    nn::Int64=lx*ly # number of spins
    nb::Int64=2*nn  # number of bonds
    mm::Int64=max(4,fld(nn,4)) #maximu cut-off
    nh::Int64=0 # number of H-operators in string 
    
    
    spin=Array{Int64}(undef,nn)
    for i=1:nn
        spin[i]=rand(-1:2:1) #initial spin configuration
    end
    bsites=Array{Int64}(undef,2,nb)
    opstring=zeros(Int64,mm)
    frstspinop=Array{Int64}(undef,nn)
    lastspinop=Array{Int64}(undef,nn)
    vertexlist=Array{Int64}(undef,4*mm)
    data1=zeros(Float64,7)
    data2=zeros(Float64,7)
    
    bsites=makelattice(nn,lx,ly,bsites) #produce bsites
    
    for i=1:isteps #for equilibration
        nh=diagonalupdate(nh,mm,opstring,spin,bsites,nb,beta)
        loopupdate(mm,nn,opstring,bsites,vertexlist,spin)
        mm=adjustcutoff(nh,mm,i,opstring,vertexlist)
    end
    
    for j=1:nbins
        
        enrg1::Float64=0.0
        enrg2::Float64=0.0
        amag2::Float64=0.0
        ususc::Float64=0.0
        for i=1:msteps
            nh=diagonalupdate(nh,mm,opstring,spin,bsites,nb,beta)
            loopupdate(mm,nn,opstring,bsites,vertexlist,spin)
            (enrg1,enrg2,amag2,ususc)=measureobservables(mm,nh,nn,lx,spin,bsites,opstring,enrg1,enrg2,amag2,ususc)
        end
        
        writeresults(nn,nb,msteps,beta,enrg1,enrg2,amag2,ususc)
    end
end

f=open("read.in","r")
    lx=parse(Int64,readline(f)) #x lattice size
    ly=parse(Int64,readline(f)) #y lattice size
    beta=parse(Int64,readline(f)) #inverse Temperature
    nbins=parse(Int64,readline(f)) # # of bins
    isteps=parse(Int64,readline(f)) #Number of MC sweeps for equilibration
    msteps=parse(Int64,readline(f)) #Number of MC sweeps in each bin
close(f)

@time mainprogram(lx,ly,beta,nbins,isteps,msteps)
