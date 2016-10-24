# Written by Marco Pingaro and Paolo Venini
#
# The script convert the figure in saved from Matlab/Octave in eps
# to figure in eps available in tex file

import os;

def CreateFile(name_input): 
    out_file = open(name_input+'.sh',"w")
    out_file.write("epstopdf %s.eps\n" % name_input)
    out_file.write("pdfcrop %s.pdf\n" % name_input)
    out_file.write("pdftops -eps %s-crop.pdf\n" % name_input)
    out_file.close()
    return;

name_input = input('Insert name of figure in eps: ') # Name of figure without extension  
CreateFile(name_input)   
os.system("sh %s.sh" % name_input)  
