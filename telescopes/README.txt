##JLINE - 26/09/2016
A directory to contain all telescope models. Each directory contains:
   - The telecope layout as read in by OSKAR
   - A template OSKAR .ini and .uvfits file for use by MWAobs_oskar_ascii.py

NOTES:
   - The array coords (east, north, height, local topocentric) defined in 
     MWA_Tools and as used by OSKAR must be negative of one another to 
     produce the same u,v,w coords. Go figure, probably something to do with
     differing hemispheres of the MWA and OSKAR??
   - The uvfits files created are tailored to work in the RTS 
     (Mitchell et al 2008). They SHOULD be industry standard but you never know
