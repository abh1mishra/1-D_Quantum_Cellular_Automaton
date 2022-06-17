using LinearAlgebra
using Luxor
using Colors
steps=15
n=5#number of qubits
pixelSize=40
canWidth=pixelSize*2^n
canHeight=pixelSize*(steps+1)

function bitarr_to_int(arr)
    return sum(arr .* (2 .^ collect(length(arr)-1:-1:0)))
end
function computeCoeffMatrix()
    coeff=zeros(Complex{Float64},2^n)
    coeff[1]=(1.0+im)
    coeff[4]=(1.0-im)
    coeff[8]=(-1.0+im)
    coeff[16]=(-1.0-im)
    coeff[32]=(-1.0-im)
    coeff/=norm(coeff)
    permMatrix=zeros(Int64,2^n,2^n)
    for i in 0:2^n-1
        bitstr=bitstring(UInt8(i))[end-n+1:end]
        bitarr=parse.(Int,split(bitstr,""))
        temp=bitarr[end]
        bitarr[end]=bitarr[2]
        bitarr[2]=temp
        permMatrix[i+1,bitarr_to_int(bitarr)+1]=1
    end
    interactionMatrix=[1.0;1.0;1.0;1.0 ;;exp(pi*0*im/4);exp(pi*im/4);exp(pi*2*im/4);exp(pi*3*im/4);;
    exp(pi*0*im/4);exp(pi*2*im/4);exp(pi*4*im/4);exp(pi*6*im/4);;exp(pi*0*im/4);exp(pi*3*im/4);exp(pi*6*im/4);exp(pi*9*im/4)]
    idMatrix=[1.0 0.;0.0 1.0]
    #generate an array with unitary acting on 1st two qubits and identity on rest
    operatorArray=[]
    push!(operatorArray,interactionMatrix)
    for i in 1:n-2
        push!(operatorArray,idMatrix)
    end
    #generate the actual QCA operators which when act will produce psi_1,psi_2...
    circOperators=[]
    for i in 1:n-1
        push!(circOperators,kron(tuple(operatorArray...)...))
        operatorArray=circshift(operatorArray,1)
    end
    stepCoeff=[]
    push!(stepCoeff,coeff)
    for i in 1:steps
        psi=[]#psi=[psi_1,psi_2...]
        for j in 1:n-1
            push!(psi,circOperators[j]*coeff)
        end
        push!(psi,permMatrix*circOperators[1]*permMatrix*coeff)
        #sum the psi_1,psi_2...
        coeff=sum(psi)/norm(sum(psi))
        push!(stepCoeff,coeff)

    end
    return stepCoeff

end
function getR(ang)
    r=0
    if 0<=ang<=120
        r=1.0
    elseif 120<=ang<=180
        r=1-((ang-120)/60)
    elseif -60<=ang<0
        r=1+((ang)/60)
    else
        r=0
    end
    return r
end
function getG(ang)
    g=0
    if (120<=ang<=180)&(-180<ang<=-120)
        g=1.0
    elseif 60<=ang<120
        g=(ang-60)/60
    elseif -120<ang<=-60
        g=1-((ang+120)/60)
    else
        g=0
    end
    return g
end
function getB(ang)
    b=0
    if -120<=ang<=0
        b=1.0
    elseif 0<ang<=60
        b=(60-ang)/60
    elseif -180<=ang<-120
        b=(180+ang)/60
    else
        b=0
    end
    return b
end

function rgbColor(c)
    ang=angle(c)*180/pi
    radius=abs(c)
    r=getR(ang)
    g=getG(ang)
    b=getB(ang)
    r=r+(1-r)*(1-radius)
    g=g+(1-g)*(1-radius)
    b=b+(1-b)*(1-radius)
    #color inversion for fun, we love black
    r=1-r
    g=1-g
    b=1-b
    color=(r,g,b)
    return color
end

function rec(ind,jnd,col=())
    sethue(RGB(col...))
    rect(O-(canWidth/2-(jnd-1)*pixelSize,canHeight/2-(ind-1)*pixelSize),pixelSize,pixelSize,:fill)
end
compMatrix=computeCoeffMatrix()
@png begin
    for i in 1:1+steps
        for j in 1:2^n
            rec(i,j,rgbColor(compMatrix[i][j]))
        end
    end
end canWidth canHeight "n=4_10000"
