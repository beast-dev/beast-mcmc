package dr.util;

import java.util.Formatter;

/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */

import java.util.Locale;

import dr.math.MathUtils;

import jebl.util.FixedBitSet;

public class AlloppMisc {
	
	
	public static String FixedBitSetasText(FixedBitSet x) {
		StringBuilder rep = new StringBuilder();
		rep.append("{");
		for (int b = 0; b < x.cardinality() + FixedBitSet.complement(x).cardinality(); ++b) {
			String comma = (b==0) ? "" : ",";
			if (x.contains(b)) {
				rep.append(comma + b);
			} else {
				rep.append(comma + " ");
			}
		}
		rep.append("}");
		return rep.toString();
	}
	
	
	
	public static String nonnegIn8Chars(double x)
	{
		StringBuilder s = new StringBuilder();
		Formatter formatter = new Formatter(s, Locale.US);
		if (x < 0) {
			formatter.format("%8s", "NA");
		} else if (x == 0.0) {
			formatter.format("%8s", "zero");
		} else if (x < 1e-3) {
			formatter.format("%8.2e", x);
		} else {
			formatter.format("%8.6f", x);
		}
		return s.toString();
	}

	
	public static String nonnegIntIn2Chars(int x)
	{
		StringBuilder s = new StringBuilder();
		Formatter formatter = new Formatter(s, Locale.US);
		if (x < 0) {
			formatter.format("%2s", "NA");
		} else {
			formatter.format("%2d", x);
		}
		return s.toString();
	}
	
	
    public static double uniformInRange(double oldx, double min, double max, double halfwidth) {
   	 assert halfwidth > 0.0;
   	 assert halfwidth < 0.5;
   	 assert min < max;
   	 double change = MathUtils.uniform(-1.0, 1.0) * (max - min) * halfwidth;
   	 double newx = oldx + change;
   	 if (newx < min) { newx = 2*min - newx; }
   	 if (newx > max) { newx = 2*max - newx; }
   	 assert newx > min;
   	 assert newx < max;
   	 return newx;
    }
	
}


