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
### script for display purposes: chemistry ####
###############################################

#load the data
t  = loadData("output/time.dat",8)
t=squeeze(t,1)
iso = loadData("output/emi_iso.dat",25)
alp = loadData("output/emi_alp.dat",25)
oh = loadData("output/OH.dat",25)
ho2 = loadData("output/HO2.dat",25)
h2so4 = loadData("output/H2SO4.dat",25)
isoprene = loadData("output/isoprene.dat",25)
alphap = loadData("output/alpha.dat",25)
elvoc = loadData("output/ELVOC.dat",25)
h= [0.0,   10,   20,   30,   40, 50,   60,   70,   80,   90, 100,  120,  140,  160,  180,  200,  230,  260,  300,  350,  400,  450,  500,  550,  600,  650,  700,  800,  900, 1000,  1100, 1200, 1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2100, 2200, 2300, 2400, 2500, 2600, 2700, 2800, 2900, 3000]


# emission rates of isoprene and alpha=pinene
ax=figure(101)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
ax[:ticklabel_format](axis="x",style="sci",format="%1.1f")
plot(t,squeeze(iso,1), linewidth=2,marker="o")
plot(t,squeeze(alp,1), linewidth=2,marker="o")
title("emission rates")
legend(["isoprene","\$\\alpha\$-pinene"])
xlabel("time [d]")
ylabel("emission rate [molec cm\$^{-3}\$ s\$^{-1}\$]")
#ax[:axes][1][:yaxis][:set_major_formatter]


# OH concentration (mapping time and height)
imshowData(102,t,h,oh,_norm=:Normalize,_vmin=minimum(oh),_vmax=maximum(oh),_edgecolors="face")
#colorbar(format="%1.1e",label="concentration [cm^{-3}]")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "OH (%1.1e,%1.1e) cm\$^{-3}\$" minimum(oh) maximum(oh)
title(s)
xlabel("time [d]")
ylabel("height [m]")


# HO2 concentration (mapping time and height)
imshowData(103,t,h,ho2,_norm=:Normalize,_vmin=minimum(ho2),_vmax=maximum(ho2),_edgecolors="face")
#colorbar(format="%1.1e",label="concentration [cm^{-3}]")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "HO2 (%1.2e,%1.2e) cm\$^{-3}\$" minimum(ho2) maximum(ho2)
title(s)
xlabel("time [d]")
ylabel("height [m]")


# H2SO4 concentration (mapping time and height)
imshowData(104,t,h,h2so4,_norm=:Normalize,_vmin=minimum(h2so4),_vmax=maximum(h2so4),_edgecolors="face")
#colorbar(format="%1.1e",label="concentration [cm^{-3}]")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "H2SO4 (%1.2e,%1.2e) cm\$^{-3}\$" minimum(h2so4) maximum(h2so4)
title(s)
xlabel("time [d]")
ylabel("height [m]")

# isoprene concentration (mapping time and height)
imshowData(105,t,h,isoprene,_norm=:Normalize,_vmin=minimum(isoprene),_vmax=maximum(isoprene),_edgecolors="face")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "isoprene (%1.2e,%1.2e) cm\$^{-3}\$" minimum(isoprene) maximum(isoprene)
title(s)
xlabel("time [d]")
ylabel("height [m]")

# alpha-pinene concentration (mapping time and height)
imshowData(106,t,h,alphap,_norm=:Normalize,_vmin=minimum(alphap),_vmax=maximum(alphap),_edgecolors="face")
#colorbar(format="%1.1e",label="concentration [cm^{-3}]")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "\$\\alpha\$-pinene (%1.2e,%1.2e) cm\$^{-3}\$" minimum(alphap) maximum(alphap)
title(s)
xlabel("time [d]")
ylabel("height [m]")

# ELVOC concentration (mapping time and height)
imshowData(107,t,h,elvoc,_norm=:Normalize,_vmin=minimum(elvoc),_vmax=maximum(elvoc),_edgecolors="face")
#colorbar(format="%1.1e",label="concentration [cm^{-3}]")
cbar=colorbar()
cbar[:set_label]("concentration [molec cm\$^{-3}\$]")
cbar[:formatter][:set_powerlimits]((-1,2))
cbar[:update_ticks]()
s = @sprintf "ELVOC (%1.2e,%1.2e) cm\$^{-3}\$" minimum(elvoc) maximum(elvoc)
title(s)
xlabel("time [d]")
ylabel("height [m]")


# OH concentration (time evolution at 10m and 50m high)
figure(108)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "OH (%1.2e,%1.2e) cm\$^{-3}\$" minimum(oh) maximum(oh)
title(s)
plot(t,squeeze(oh[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(oh[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("oh [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])

# HO2 concentration (time evolution at 10m and 50m high)
figure(109)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "HO2 (%1.2e,%1.2e) cm\$^{-3}\$" minimum(ho2) maximum(ho2)
title(s)
plot(t,squeeze(ho2[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(ho2[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("oh [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])

# H2SO4 concentration (time evolution at 10m and 50m high)
figure(110)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "H2SO4 (%1.2e,%1.2e) cm\$^{-3}\$" minimum(h2so4) maximum(h2so4)
title(s)
plot(t,squeeze(h2so4[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(h2so4[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("oh [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])

# isoprene concentration (time evolution at 10m and 50m high)
figure(111)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "isoprene (%1.2e,%1.2e) cm\$^{-3}\$" minimum(isoprene) maximum(isoprene)
title(s)
plot(t,squeeze(isoprene[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(isoprene[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("isoprene [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])

# alpha-pinene concentration (time evolution at 10m and 50m high)
figure(112)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "\$\\alpha\$-pinene (%1.2e,%1.2e) cm\$^{-3}\$" minimum(alphap) maximum(alphap)
title(s)
plot(t,squeeze(alphap[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(alphap[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("\$\\alpha\$-pinene [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])

# ELVOC concentration (time evolution at 10m and 50m high)
figure(113)
ax=axes()
ax[:ticklabel_format](axis="y",style="sci",scilimits=(2,2))
s = @sprintf "ELVOC (%1.2e,%1.2e) cm\$^{-3}\$" minimum(elvoc) maximum(elvoc)
title(s)
plot(t,squeeze(elvoc[2,:],1), linewidth=2,marker="o")
plot(t,squeeze(elvoc[6,:],1), linewidth=2,marker="o")
xlabel("time [d]")
ylabel("ELVOC [molec cm\$^{-3}\$]")
legend(["10 m","50 m"])
