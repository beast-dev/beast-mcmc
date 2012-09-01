package dr.util;

import java.util.Formatter;
import java.util.Iterator;

/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */

import java.util.Locale;

import dr.evolution.tree.NodeRef;
import dr.evolution.tree.SimpleTree;
import dr.evolution.util.Taxon;
import dr.math.MathUtils;

import jebl.util.FixedBitSet;

public class AlloppMisc {
	
	
	public static String FixedBitSetasText(FixedBitSet x) {
        if (x == null) {
            return "{?}";
        }
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
		} else if (x < 9.999) {
			formatter.format("%8.5f", x);
		} else if (x < 99.99) {
            formatter.format("%8.4f", x);
        } else if (x < 999.9) {
            formatter.format("%8.3f", x);
        } else if (x < 9999) {
            formatter.format("%8.2f", x);
        } else {
            formatter.format("%8.0f", x);
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
	

	public static String SimpleNodeAsText(SimpleTree stree, NodeRef node) {
		String s = "" + node.getNumber() + " ";
		while (s.length() < 3) { s += " "; }
		int nch = stree.getChildCount(node);
		if (nch> 0) {
			assert(nch==2);
			s += stree.getChild(node, 0).getNumber();
			while (s.length() < 6) { s += " "; }
			s += stree.getChild(node, 1).getNumber();
		}
		while (s.length() < 9) { s += " "; }
		Taxon tx = stree.getNodeTaxon(node);
		String taxonid = "*";
		if (tx != null) {
			taxonid = tx.getId();
			if (taxonid == null || taxonid.length() == 0) {
				taxonid = "*";
			}

		}
		s += taxonid;
		while (s.length() < 20) { s += " "; }

		stree.getNodeHeight(node);
		s += " height=";
		s += stree.getNodeHeight(node);
		Iterator iter = stree.getNodeAttributeNames(node);
		if (iter != null) {
			while (iter.hasNext()) {
				String name = (String) iter.next();
				s += " ";
				s += name;
				s += "=";
				s += stree.getNodeAttribute(node, name);
			}
		}
		return s;
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


