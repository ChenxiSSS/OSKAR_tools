from numpy import log,pi,arcsin,sin,cos,isnan,nan,sqrt
from optparse import OptionParser
from sys import exit

parser = OptionParser()

parser.add_option('-f','--extrap_freq', help='Frequency to extrap sources to (in Hz)')
parser.add_option('-s','--srclist', help='Name of srclist to convert')
parser.add_option('-o','--osmname', help='Name for output sky model name (include .osm)')

options, args = parser.parse_args()

extrap_freq = float(options.extrap_freq)

##Convert catalogue/OSKAR gaussians into RTS gaussians
fudge = sqrt(pi**2 / (2*log(2)))

##Class to store source information with - set to lists
##to store component info in the same place
class rts_source():
    def __init__(self):
        self.name = ''
        self.ras = []
        self.decs = []
        self.freqs = []
        self.fluxs = []
        self.SIs = []
        self.extrap_fluxs = []
        self.pas = []
        self.majors = []
        self.minors = []
        
def find_SI(freqs,fluxs,extrap_freq):
    '''f1/f2 = (nu1/n2)**alpha
       alpha = ln(f1/f2) / ln(nu1/nu2)
       f1 = f2*(nu1/nu2)**alpha'''
    alpha = log(fluxs[0]/fluxs[1]) / log(freqs[0]/freqs[1])
    extrap_flux = fluxs[0]*(extrap_freq/freqs[0])**alpha
    return alpha,extrap_flux

def create_sources(split_source):
    
    source = rts_source()
    
    ##Find the primary source info - even if no comps, this will isolate
    ##primary source infomation
    primary_info = split_source.split('COMPONENT')[0].split('\n')
    primary_info = [info for info in primary_info if info!='']
    meh,prim_name,prim_ra,prim_dec = primary_info[0].split()
    
    ##Put in to the source class
    source.name = prim_name
    source.ras.append(float(prim_ra)*15.0)
    source.decs.append(float(prim_dec))
    ##Find the fluxes and append to the source class
    prim_freqs = []
    prim_fluxs = []
    yes_gauss = False
    for line in primary_info:
        if 'FREQ' in line:
            prim_freqs.append(float(line.split()[1]))
            prim_fluxs.append(float(line.split()[2]))
        elif 'GAUSSIAN' in line:
            yes_gauss = True
            meh,pa,major,minor = line.split()
            source.pas.append(float(pa))
            ##Change from arcmin to arcsec (*60)
            ##Convert into OSKAR major ( / fudge)
            source.majors.append((float(major)*60.0) / fudge)
            source.minors.append((float(minor)*60.0) / fudge)
            
    if yes_gauss:
        pass
    else:
        source.pas.append(nan)
        source.majors.append(nan)
        source.minors.append(nan)
            
    source.freqs.append(prim_freqs)
    source.fluxs.append(prim_fluxs)
    
    ##Split all info into lines and get rid of blank entries
    lines = split_source.split('\n')
    lines = [line for line in lines if line!='']
    ##If there are components to the source, see where the components start and end
    comp_starts = [lines.index(line) for line in lines if 'COMPONENT' in line and 'END' not in line]
    comp_ends = [i for i in xrange(len(lines)) if lines[i]=='ENDCOMPONENT']
    
    ##For each component, go through and find ra,dec,freqs and fluxs
    for start,end in zip(comp_starts,comp_ends):
        freqs = []
        fluxs = []
        yes_gauss = False
        for line in lines[start:end]:
            if 'COMPONENT' in line:
                source.ras.append(float(line.split()[1]))
                source.decs.append(float(line.split()[2]))
            elif 'FREQ' in line:
                freqs.append(float(line.split()[1]))
                fluxs.append(float(line.split()[2]))
            elif 'GAUSSIAN' in line:
                yes_gauss = True
                meh,pa,major,minor = line.split()
                source.pas.append(float(pa))
                ##Change from arcmin to arcsec (*60)
                ##Convert into OSKAR major ( / fudge)
                source.majors.append((float(major)*60.0) / fudge)
                source.minors.append((float(minor)*60.0) / fudge)
                
        if yes_gauss:
            pass
        else:
            source.pas.append(nan)
            source.majors.append(nan)
            source.minors.append(nan)
                
        source.fluxs.append(fluxs)
        source.freqs.append(freqs)
    
    ##Mimic the way that the RTS would try to extrapolate the flux of the source
    ##to the observational frequency
    for freqs,fluxs in zip(source.freqs,source.fluxs):
        print freqs,fluxs
        ##If one flux, assume -0.7 like the RTS
        if len(freqs)==1:
            SI = -0.7
            ext_flux = fluxs[0]*(extrap_freq/freqs[0])**SI
        ##If extrapolating below known freqs, choose two lowest frequencies
        elif min(freqs)>extrap_freq:
            SI,ext_flux = find_SI([freqs[0],freqs[1]],[fluxs[0],fluxs[1]],extrap_freq)
        ##If extrapolating above known freqs, choose two highest frequencies
        elif max(freqs)<extrap_freq:
            SI,ext_flux = find_SI([freqs[-2],freqs[-1]],[fluxs[-2],fluxs[-1]],extrap_freq)
        ##Otherwise, choose the two frequencies above and below, and extrap between them
        else:
            for i in xrange(len(freqs)-1):
                if freqs[i]<=extrap_freq and freqs[i+1]>extrap_freq or freqs[i]<extrap_freq and freqs[i+1]>=extrap_freq:
                    
                    #print source.name, freqs[i],freqs[i+1],fluxs[i],fluxs[i+1]
                    
                    SI,ext_flux = find_SI([freqs[i],freqs[i+1]],[fluxs[i],fluxs[i+1]],extrap_freq)
        
        source.SIs.append(SI)
        source.extrap_fluxs.append(ext_flux)
        
        if isnan(ext_flux):
            print "Extrapolated flux is NAN, you gonna have problems",exit(0)
        
    #sources.append(source)
    return source
        
rts_srcs = open(options.srclist,'r').read().split('ENDSOURCE')
del rts_srcs[-1]

#print rts_srcs

oskar_outfile = open(options.osmname,'w+')
#sources = []
##Go through all sources in the source list, gather their information, extrapolate
##the flux to the central frequency and weight by the beam at that position
for split_source in rts_srcs:

    rts_src = create_sources(split_source)
    ##Format of oskar osm file - put in the extrapolated flux and SI at the observational frequency
    #  RA,    Dec,   I,    Q,    U,    V,   freq0, spix,  RM,      maj,      min,      pa
    # (deg), (deg), (Jy), (Jy), (Jy), (Jy), (Hz), (-), (rad/m^2), (arcsec), (arcsec), (deg)
    
    for i in xrange(len(rts_src.ras)):
        if isnan(rts_src.pas[i]) == True:
            oskar_outfile.write('%.10f %.10f %.10f 0 0 0 %.5f %.3f 0.0 0 0 0\n' %(rts_src.ras[i],rts_src.decs[i],rts_src.extrap_fluxs[i],extrap_freq,rts_src.SIs[i]))
        else:
            oskar_outfile.write('%.10f %.10f %.10f 0 0 0 %.5f %.3f 0.0 %.2f %.2f %.1f\n' %(rts_src.ras[i],rts_src.decs[i],rts_src.extrap_fluxs[i],extrap_freq,rts_src.SIs[i],rts_src.majors[i],rts_src.minors[i],rts_src.pas[i]))
    
oskar_outfile.close()    
    

    
