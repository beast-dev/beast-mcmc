from __future__ import division
from math import exp
from scipy import array
from scipy.stats import scoreatpercentile

__all__ = ["plotFromCSV", "plotFromAll"]

import pylab

def drawBetween(x, yl, yh, col, lw, alpha = 1, plot = pylab.plot) :
  fx = pylab.concatenate( (x,x[::-1]) )
  fy = pylab.concatenate( (yh,yl[::-1]) )

  # probably does not work with log??
  p = pylab.fill(fx, fy, facecolor=col, lw = 0, alpha = alpha)
  if lw :
    plot(x, yl, x, yh, aa = 1, alpha = alpha, lw = lw, color='k')

def locateIndex(possibleColNames, columns, prefixOK = False, noCase = False) :

  if noCase :
    columns = [c.lower() for c in columns]
    
  for s in possibleColNames :
    for k,c in enumerate(columns):
      if prefixOK:
        if c.startswith(s) :
          return k
      else:
        if s == c:
          return k

  return -1

def plotFromCSV(fName, trued, truecol = 'b', medcol = 'k', 
                pntcol = None, hpdcol = None, alpha = 1, ls = "--",
                logy = False, doHPD = True, doHist = False, xScale = 1, yScale = 1, mainlw = 2,
                xLabel = "Time", yLabel = "Population", medLabel = "median", tint = None,
                hpdOutline = 1, labelProps = dict()) :
  l = file(fName).read()
  lines = l.split('\n')

  while True:
    if lines[0][:4].lower() == "time" :
      break
    del lines[0]

  h = lines[0]
  plot = pylab.semilogy if logy else pylab.plot
  
  for n,x in enumerate(lines[1:]) :
    if len(x) == 0 or x[0].isalpha() :
      break

  sep = ','
  if sep not in h:
    sep = '\t'
    if sep not in h:
      raise RuntimeError("can't determine seperator")
    
  h = h.split(sep)

  it = locateIndex(("time",),   h, noCase=True)
  im = locateIndex(("median",), h, noCase=True)
  ib = locateIndex(("bins",), h, noCase=True) if doHist else -1
  
  if it < 0 or im < 0:
    raise RuntimeError("can't determine column positions" + str(h))

  if doHist and ib < 0 :
    raise RuntimeError("can't determine column position of bins")
    
  if doHPD :
    iHPD = array(range(len(h))).compress(
      [x.startswith("hpd") or x.startswith("cpd") for x in h])
    
    if len(iHPD) == 0 :
      hpdlow = locateIndex(("Lower",), h)
      hpdhigh = locateIndex(("Upper",), h)
      iHPD = array((hpdlow, hpdhigh))
      #print iHPD
  else :
    iHPD = []
    
  g = [x.split(sep) for x in lines[1:n+1]]

  if ib > 0 :
    binsData = [int(a[ib]) for a in g if a[ib] != ""]
    
  lastTime = len(g)
  while g[lastTime-1][0] == "" :
    lastTime -= 1
  g = g[:lastTime]
  
  x = [float(a[it]) * xScale for a in g]
  m = [float(a[im]) * yScale for a in g]

  hpd = list()
  for i in iHPD:
    hpd.append([float(a[i]) * yScale for a in g])
  #pl = [a[ipl] for a in g]

  hpd2valfunc = lambda x : 1 - (pylab.exp(x) - 1)/(exp(1)-1)
  for i in range(0, len(iHPD), 2) :
    col = h[iHPD[i]]
    if col.startswith("hpd lower") or col.startswith("cpd lower"):
      s = h[iHPD[i]][len("hpd lower"):]
      p = 1 - float(s if len(s) else 95)/100.
      c = hpd2valfunc(p)
      c = [c,c,c]
      if tint is not None :
        for z in tint:
          c[z] = 0
        
      #print i,c, s,p, h[iHPD[i]]
    else :
      assert col == "Lower"
      c = hpd2valfunc(.95)
      c = (c,c,c)

    if hpdcol :
      c = [a1*a2 for a1,a2 in zip(c, hpdcol)]

    drawBetween(x, hpd[i], hpd[i+1], c, hpdOutline, alpha = alpha, plot = plot)

  xl = x[-1]

  if ls :
    plot(x, m, "", linewidth = mainlw, ls = ls, label = medLabel, color = medcol)
  else :
    plot(x, m, "", linewidth = mainlw, label = medLabel, color = medcol)
    
  if pntcol:
    plot(x, m, pntcol + 'o')

  if ib > 0 :
    y0,y1 = pylab.ylim()
    if logy :
      mx = scoreatpercentile(binsData, 95)
      z1 = exp(0.8 * pylab.log(y0) + 0.2 * pylab.log(y1))
      # allow one more for histogram
      y0 /= 10.
      v = [max(z1 * (pylab.exp(z/mx) - 1)/(pylab.exp(1)-1),y0) for z in binsData]
    else :
      mn, mx = min(binsData), scoreatpercentile(binsData, 95)
      sc = (y1 - y0) * 0.1
      v = [y0 + sc*(z-mn)/(mx-mn) for z in binsData]
      print sc
    #print v
    ww = x[-1]/len(v)
    pylab.bar(pylab.linspace(0, x[-1], len(v)), v, width = ww, alpha = .2, color = "orange")
    pylab.ylim( (y0, y1) )
    
  if xLabel:
    pylab.xlabel(xLabel, **labelProps)
  if yLabel:
    pylab.ylabel(yLabel, **labelProps)  
  pylab.xlim( (0, xl) )


def plotFromAll(fName,
                #trued, truecol = 'b',
                medcol = 'g', 
                #pntcol = None, hpdcol = None,
                alpha = .1,
                # ls = "--",
                logy = False,
                #doHPD = True, doHist = False,
                xScale = 1, yScale = 1,
                mainlw = 1,
                xLabel = "Time", yLabel = "Population",
                # medLabel = "median",
                #tint = None,
                #hpdOutline = 1,
                labelProps = dict()) :
  fall = file(fName)
  xs = [xScale * float(x) for x in fall.readline().strip().split()]

  plot = pylab.semilogy if logy else pylab.plot

  for l in fall:
    if yScale != 1.0 :
      ys = [yScale * float(x) for x in l]
    else :
      ys = [float(x) for x in l.strip().split()]
      
    plot(xs, ys, "", linewidth = mainlw, color = medcol, alpha = alpha)
    
  if xLabel:
    pylab.xlabel(xLabel, **labelProps)
  if yLabel:
    pylab.ylabel(yLabel, **labelProps)  
  pylab.xlim( (0, xs[-1]) )

 # pylab.legend(loc='best')
