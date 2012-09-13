#!/usr/bin/env python

import sys, os.path, math, fnmatch
from glob import glob

import optparse


from popGraphUtil import plotFromCSV, plotFromAll 

parser = optparse.OptionParser(" [options] csv-file chart-file")

parser.add_option("", "--xlim", dest="xlim", help="cut off X-axis at this point", default = None)
parser.add_option("", "--ylim", dest="ylim", help="cut off Y-axis at this point", default = None)

parser.add_option("", "--logy", dest="logy", action="store_true",
                  help="Log scale for Y axis", default = False)
parser.add_option("", "--yscale", dest="yscale", help="Y-axis scale factor", default = 1)

parser.add_option("", "--width", dest="width",
                  help="figure width. Integral value with units: 50mm 2cm 3 (inches)", default = None)
# parser.add_option("", "--ms", dest="ms", help="", default = None)
parser.add_option("", "--lw", dest="lw", help="Line width", default = None)

parser.add_option("", "--font", dest="font", help="name of font for figure text ", default = None)
parser.add_option("", "--fontsize", dest="fontsize", help="font size of figure text", default = None)

# parser.add_option("", "--axes", dest="axesSize", help="", default = None)
parser.add_option("", "--ticks", dest="ticklabelsize",
                  help="font size of ticks labels ", default = None)

parser.add_option("", "--nxticks", dest="nxticks",
                  help="number of X-axis ticks", default = None)

parser.add_option("", "--title", dest="title",
                  help="Figure title", default = None)

parser.add_option("", "--hist", dest="hist", action="store_true",help="", default = False)

parser.add_option("", "--alldemo", dest="alldfile",
                  help="plot all demographic functions in this file",
                  default = None)

parser.add_option("-a", "--alphaout", dest="alpha", help="transparancy value of outline.", default = 1)


parser.add_option("", "--alpha", dest="alldalpha",
                  help="transparancy value to use when plotting all" + 
                  " demographic. 1 - no transparancy, 0 fully transparent.", default = 0.1)

parser.add_option("", "--ratio", dest="ratio",
                  help="height/width ratio of figure.", default = 0.75)

options, args = parser.parse_args()
if len(args) != 2 :
  print >> sys.stderr, "usage:", sys.argv[0], "csv-file", "chart-file"
  sys.exit(1)
  
name = args[0]


trueDemo = None

plotOptionsDict = { 'alpha' : float(options.alpha),
                    'logy' : options.logy,
                    'doHist': options.hist }
if options.lw :
  plotOptionsDict['mainlw'] = float(options.lw)
  plotOptionsDict['hpdOutline'] = float(options.lw)/2

labelsFont = None
if options.font :
  import matplotlib.font_manager
  labelsFont = matplotlib.font_manager.FontProperties(options.font)
  if labelsFont.get_name() != options.font :
    print >> sys.stderr, "*warning:", labelsFont.get_name(),"!=",options.font
  if options.fontsize :
    labelsFont.set_size(float(options.fontsize))

import pylab

def convertToInches(w) :
  if w[-2:] == 'mm' :
    return int(w[:-2]) / 25.4
  if w[-2:] == 'cm' :
    return int(w[:-2]) / 2.54
  return int(w)

if options.width is None :
  fig = pylab.figure()
else :
  w = convertToInches(options.width)
  h = w * float(options.ratio)
  fig = pylab.figure(figsize=(w,h))

if labelsFont :
  labelFontDict = {'fontproperties': labelsFont}
  plotOptionsDict['labelProps'] = labelFontDict

if options.alldfile:
  pylab.ioff()
  plotFromAll(options.alldfile, yScale = float(options.yscale),
              logy = options.logy, alpha = float(options.alldalpha))

plotFromCSV(name, trueDemo, yScale = float(options.yscale), **plotOptionsDict)
  
if options.xlim :
  pylab.xlim((0, float(options.xlim)))

if options.ylim :
  pylab.ylim((0, float(options.ylim)))

if options.title :
  pylab.title(options.title)
  
pylab.legend(loc='best')

if options.nxticks :
  from matplotlib.ticker import MaxNLocator
  pylab.gca().xaxis.set_major_locator(MaxNLocator(int(options.nxticks)))

if labelsFont :
  ltext = pylab.gca().get_legend().get_texts()
  for l in ltext :
    pylab.setp(l, fontproperties = labelsFont)

if options.ticklabelsize :
  s = float(options.ticklabelsize)
  if labelsFont :
    fp = matplotlib.font_manager.FontProperties(labelsFont.get_name())
    fp.set_size(s)
    fp = {'fontproperties' : fp}
  else :
    fp = dict()

  for p in ('xticklabels', 'yticklabels') :
    l = pylab.getp(pylab.gca(), p)
    pylab.setp(l, fontsize=s, **fp)

if options.alldfile:
  pylab.ion()
  
pylab.savefig(args[1], dpi=300)

