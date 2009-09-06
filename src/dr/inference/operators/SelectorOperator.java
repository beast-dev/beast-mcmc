package dr.inference.operators;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.util.List;
import java.util.ArrayList;

/**
 * The code is much more elegant in Python. I definitly don't have to write my own max on an list of integers (yikes).
 *
 * @author Joseph Heled
 *         Date: 4/09/2009
 */
public class SelectorOperator extends SimpleMCMCOperator {
    public static String SELECTOR_OPERATOR = "selectorOperator";
    private final Parameter selector;

    private final int[] np;

    SelectorOperator(Parameter selector) {
        this.selector = selector;
        final int len = selector.getSize();
        np = new int[len +1];

        for(int l = 0; l < np.length; ++l) {
            np[l] = npos(len, l);
        }
    }
    
    public String getOperatorName() {
        return SELECTOR_OPERATOR + "(" + selector.getParameterName() + ")";
    }

    public double doOperation() throws OperatorFailedException {

        final int[] s = vals();
        final List<Integer> poss = movesFrom(s);
        final int i = MathUtils.nextInt(poss.size()/2);

        final int[] y = new int[s.length];
        System.arraycopy(s, 0, y, 0, s.length);
        final Integer p = poss.get(2 * i);
        y[p] = poss.get(2*i+1);

        double hr = count1sr(s, y);
        hr *= (double)(poss.size()*np[max(s)])/(movesFrom(y).size() * np[max(y)]);

        selector.setParameterValue(p, y[p]);
        
        return Math.log(hr);
    }

    public String getPerformanceSuggestion() {
        return null;
    }

    private int[] vals() {
        int[] v = new int[selector.getSize()];
        for(int k = 0; k < v.length; ++k) {
            v[k] = (int)(selector.getParameterValue(k) +0.5);
        }
        return v;
    }

    private int npos(int s, int m) {
        return npos(s, m, 1);
    }

    private int npos(int s, int m, int mn) {
        if( m == 0 || s == 0 ) {
            return 1;
        }

        int tot = 0;
        for(int k = mn; k < 1+s/m; ++k) {
            final int r = s - k*(m+1);
            if(r < 0 ) {
                break;
            }
            tot += npos(r, m-1, 0);
        }
        return tot;
    }

    private int sum(int[] s) {
        int sum = 0;
        for(int si : s) {
            sum += si;
        }
        return sum;
    }

    private int max(int[] s) {
        int mx = s[0];

        for(int k = 1; k < s.length; ++k) {
            if( mx < s[k] ) {
                mx = s[k];
            }
        }
        return mx;
    }

    private int[] counts(int[] s, int mx) {
        int[] c = new int[mx+1];
        for(int si : s) {
            c[si]++;
        }
        return c;
    }

    private List<Integer> movesFrom(int[] s) {
        final int mx = max(s);

        final int[] counts = counts(s, mx);
        final List<Integer> opt = new ArrayList<Integer>(5);

        for(int k = 0; k < s.length; ++k) {
            final int si = s[k];
            if( si<mx && ((counts[si] == 1) || counts[si] == counts[si+1]) ) {
                // only or breaks order -> no moves
            } else {
                for(int x = 0; x < mx+1; ++x) {
                    if(x == si) {
                        continue;
                    }
                    if( (x > si && counts[si] - 1 >= counts[x] + 1 && counts[x-1] >= counts[x]+1)
                            ||
                            (x < si && (x > 0 && counts[x]+1 <= counts[x-1] || x == 0)) ) {
                        opt.add(k);
                        opt.add(x);
                    }
                }
                if( counts[si] > 1) {
                    opt.add(k);
                    opt.add(mx+1);
                }
            }
        }
        return opt;
    }

  private long choose(int n, int k) {
      double r = 1;
      while( n > k ) {
          r *= n;
          r /= (n-k);
          --n;
      }
      return (long)(r+0.5);
  }

   private long[] count1l(int[] ls) {
       int l = sum(ls);
       int i = 0;
       long[] r = new long[ls.length];

       while(l > 0) {
           r[i] = choose(l, ls[i]);

           l -= ls[i];
           i += 1;
       }
       return r;
   }

    private double count1sr(int[] x, int[] y) {
        long[] r1 = count1l(counts(x, max(x)));
        long[] r2 = count1l(counts(y, max(y)));

        int k = Math.min(r1.length, r2.length);
        double r = 1;
        for(int i=0; i < k; ++i) {
            r *= r1[i];
            r /= r2[i];
        }
        for(int i=k; i < r1.length; ++i) {
            r *= r1[i];
        }
        for(int i=k; i < r2.length; ++i) {
            r /= r2[i];
        }
        return r;
    }


    public static dr.xml.XMLObjectParser PARSER = new dr.xml.AbstractXMLObjectParser() {

        public String getParserName() {
            return SELECTOR_OPERATOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            final Parameter parameter = (Parameter) xo.getChild(Parameter.class);
            return new SelectorOperator(parameter);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return ".";
        }

        public Class getReturnType() {
            return SelectorOperator.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[] {
                    new ElementRule(Parameter.class),
            };
        }
    };
}
