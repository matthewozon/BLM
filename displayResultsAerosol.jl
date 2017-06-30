    using PyPlot
# define some useful functions

# conversion of bulk string to double
function convertToFloat64(str::Array{ByteString,1},len::Int64)
    # str: string data
    # len: length of each datum in character
    
    #compute sizes
    numTimeStep=round(Int64,size(str,1))
    numLevel=round(Int64,length(str[1])/len)
    #allocate a 2D array for the data
    strDouble=Array(Cdouble,numLevel,numTimeStep)

    #convert each element
    for i in 1:numTimeStep
        for j in 1:numLevel
            idxS=1+len*(j-1)
            idxE=len*j
            strDouble[j,i]=parse(Cdouble,str[i][idxS:idxE])
        end
    end

    #return
    strDouble
end

# load data from file and return a 2D array of double
function loadData(fileName::ASCIIString,len::Int64)
    # fileName: file name containing the organized data
    # len:      length of each datum in character

    # load the data in a string array
    f=open(fileName,"r") # should check if the file is open (isopen)
    dataSTR=readlines(f)           # should check if there actually are data (isempty)
    close(f)

    # convert the string data to Float64 array
    convertToFloat64(dataSTR,len)
end

# display data on the actual grid (taking into account the pixel coordinates)
function imshowData(figNum::Int64,_t::Array{Cdouble,1},_h::Array{Cdouble,1},Z::Array{Cdouble,2};_norm=:NoNorm,_vmin=0.0,_vmax=1.0,_edgecolors="face",_shading="None")
    # get the focus on the figure or create a figure
    figure(figNum)

    # get the current axis
    ax = axes()

    # display image if possible
    if (length(_t),length(_h))==size(Z) || (length(_h),length(_t))==size(Z)
        pcolormesh(repmat(_t,1,length(_h))',repmat(_h,1,length(_t)),Z,edgecolors=_edgecolors,shading=_shading,norm=matplotlib[:colors][:Normalize](vmin=_vmin,vmax=_vmax)) #shading="gouraud"
    else
        @sprintf "Cannot display the image because of a size mismatch: (%i,%i)!=%i,%i" size(Z,1) size(Z,2) length(_t) length(_h)
    end

    #return axis handler
    ax
end


###############################################
#### script for display purposes: aerosol #####
###############################################

#load the data
t        = loadData("output/time.dat",8)
diameter = loadData("output/diameters.dat",25)
conc     = loadData("output/particle.dat",25)
pn       = loadData("output/PN.dat",25)
pm       = loadData("output/PM.dat",25)
cs       = loadData("output/CS.dat",25)
hno3     = loadData("output/HNO3.dat",25)
h        = loadData("output/h.dat",25)

# squeeze the 1D data
t = squeeze(t,1)
h = squeeze(h,2)
diameter=squeeze(diameter,2)
volume=(pi/6.0)*diameter.^3
#h= [0.0,   10,   20,   30,   40, 50,   60,   70,   80,   90, 100,  120,  140,  160,  180,  200,  230,  260,  300,  350,  400,  450,  500,  550,  600,  650,  700,  800,  900, 1000,  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]


# size distribution at the beginning and at the end of the simulation
figure(201)
loglog(10.0^9*diameter[1:end-1],10.0^(-6)*conc[1:end-1,1])
loglog(10.0^9*diameter[1:end-1],10.0^(-6)*conc[1:end-1,end])
ylim(1.0,1000000.0)
grid(true,which="major",ls="-")
grid(true,which="minor",ls="-",alpha=0.5)
title("size distribution")
xlabel("diameter [nm]")
ylabel("concentration [cm\$^{-3}\$]")
legend(["t=0d","t=5d"])


# condensation sink of H2SO4 wrt time
figure(202)
plot(t,squeeze(cs[2,:],1))
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(-1,2))
xlabel("time [d]")
ylabel("cs [Hz]")
title("condensation sink: first layer")

# particle number concentration  wrt time
figure(203)
plot(t,squeeze(pn[2,:],1)*10.0^(-6))
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(-1,2))
xlabel("time [d]")
ylabel("pn [cm\$^{-3}\$]")
title("particle number concentration: first layer")

# particle mass concentration wrt time
figure(204)
plot(t,squeeze(pm[2,:],1)*10.0^9)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(-1,2))
xlabel("time [d]")
ylabel("pm [\$\\mu\$g cm\$^{-3}\$]")
title("particle mass concentration: first layer")

# size distribution
figure(205)
plot(pn[:,end]*10.0^(-6),h)
ax=axes()
ax[:ticklabel_format](axis="both",style="sci",scilimits=(-1,2))
#ax[:ticklabel_format](axis="x",style="sci",scilimits=(-1,2))
title("particle number: last day")
ylabel("height [m]")
xlabel("pn [cm\$^{-3}\$]")

# mass distribution
figure(206)
plot(pm[:,end]*10.0^9,h)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(-1,2))
ylabel("height [m]")
xlabel("pm [\$\\mu\$g cm\$^{-3}\$]")
title("particle mass: last day")

# size distribution 
imshowData(207,t,10.0^9*diameter,log10(10.0^(-6)*conc),_norm=:Normalize,_vmin=0.0,_vmax=3.18,_edgecolors="face")
yscale("log")
cbar=colorbar()
cbar[:set_label]("concentration [cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
mini=minimum(10.0^(-6)*conc)
maxi=maximum(10.0^(-6)*conc)
s = @sprintf "size distribution (%1.2e,%1.2e) cm\$^{-3}\$" mini maxi
title(s)
xlabel("time [d]")
ylabel("height [m]")

# total particle number concentration (time and height mapping)
mini=minimum(10.0^(-6)*pn)
maxi=maximum(10.0^(-6)*pn)
imshowData(208,t,h,10.0^(-6)*pn,_norm=:Normalize,_vmin=mini,_vmax=maxi,_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("concentration [cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "particle number (%1.2e,%1.2e) cm\$^{-3}\$" mini maxi
title(s)
xlabel("time [d]")
ylabel("height [m]")

# total particle number concentration (time and height mapping)
mini=minimum(10.0^9*pm)
maxi=maximum(10.0^9*pm)
imshowData(209,t,h,10.0^9*pm,_norm=:Normalize,_vmin=mini,_vmax=maxi,_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("mass concentration [g cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "particle mass (%1.2e,%1.2e) \$\\mu\$g cm\$^{-3}\$" mini maxi
title(s)
xlabel("time [d]")
ylabel("height [m]")

# H2SO4 condensation sink (time and height mapping)
mini=minimum(cs)
maxi=maximum(cs)
imshowData(210,t,h,cs,_norm=:Normalize,_vmin=mini,_vmax=maxi,_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("frequency [Hz]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "condensation sink (%1.2e,%1.2e) Hz" mini maxi
title(s)
xlabel("time [d]")
ylabel("height [m]")



# HNO3 concentration (mapping time and height)
imshowData(211,t,h,hno3,_norm=:Normalize,_vmin=minimum(hno3),_vmax=maximum(hno3),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "HNO3 (%1.2e,%1.2e) cm\$^{-3}\$" minimum(oh) maximum(oh)
title(s)
xlabel("time [d]")
ylabel("height [m]")
