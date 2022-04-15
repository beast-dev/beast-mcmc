/*
 * SkewTDistribution.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.inference.distribution;

import dr.math.UnivariateFunction;
import dr.math.distributions.Distribution;
import dr.util.Author;
import dr.util.Citable;
import dr.util.Citation;
import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.AbstractContinuousDistribution;

import java.util.Collections;
import java.util.List;

/**
 * The Skew-T distribution in Wilkinson et al, 2011
 *
 * @author Xiang Ji
 * @author Marc Suchard
 */
public class SkewTDistribution extends AbstractContinuousDistribution implements Distribution, Citable {

    private final double location;
    private final double scale;
    private final double shape;
    private final double df;

    public SkewTDistribution (double location,
                              double scale,
                              double shape,
                              double df) {
        this.location = location;
        this.scale = scale;
        this.shape = shape;
        this.df = df;
    }

    @Override
    public double pdf(double x) {
        // PDFSkewT function in paml 4.8
        final double z = (x-location)/scale;
        final double lnghalf=0.57236494292470008707;    /* log{G(1/2)} = log{sqrt(Pi)} */

        final double lngv = LnGamma(df/2);
        final double lngv1 = LnGamma((df+1)/2);
        final double lnConst_pdft = lngv1 - lngv - 0.5*Math.log(Math.PI*df);
        final double lnbeta_cdft = lngv1 + lnghalf - lngv - Math.log(df/2);  /* log{ B((df+1)/2, 1/2) }  */

        final double pdf = 2/scale * PDFt(z, 0, 1, df, lnConst_pdft)
                * CDFt(shape*z*Math.sqrt((df+1)/(df+z*z)), 0, 1, df+1, lnbeta_cdft);
        return pdf;
    }

    private double LnGamma(double x) {
       /* LnGamma function in paml 4.8
       returns ln(gamma(x)) for x>0, accurate to 10 decimal places.
       Stirling's formula is used for the central polynomial part of the procedure.

       Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
       Communications of the Association for Computing Machinery, 9:684
        */
        assert x > 0;
        double f=0, fneg=0, z, lng;
        int nx=(int)x;

        if((double)nx==x && nx>=0 && nx<=11)
            lng = Math.log((double)factorial(nx-1));
        else {
            if(x<=0) {
                throw new RuntimeException("Undefined negative value for Gamma function.");
            }
            if (x<7) {
                f = 1;
                z = x-1;
                while (++z<7)
                    f *= z;
                x = z;
                f = -Math.log(f);
            }
            z = 1/(x*x);
            lng = fneg + f + (x-0.5)*Math.log(x) - x + 0.918938533204673
                    + (((-0.000595238095238*z + 0.000793650793651)*z - 0.002777777777778)*z + 0.083333333333333)/x;
        }
        return  lng;
    }

    private double factorial(int n) {
        long f=1, i;
        if (n>11) throw new RuntimeException("n>10 in factorial");
        for (i=2; i<=(long)n; i++) f *= i;
        return (f);
    }

    private double PDFt (double x, double loc, double scale, double df, double lnConst)
    {
        /* PDFt function from paml 4.8
        PDF of t distribution with lococation, scale, and degree of freedom
         */
        double z = (x-loc)/scale, lnpdf=lnConst;

        if(lnpdf==0) {
            lnpdf = LnGamma((df+1)/2) - LnGamma(df/2) - 0.5*Math.log(Math.PI*df);
        }
        lnpdf -= (df+1)/2 * Math.log(1+z*z/df);
        return Math.exp(lnpdf)/scale;
    }

    private double CDFt (double x, double loc, double scale, double df, double lnbeta)
    {
        /* CDFt function from paml 4.8
        CDF of t distribution with location, scale, and degree of freedom
         */
        double z = (x-loc)/scale, cdf;
        double lnghalf = 0.57236494292470008707;  /* log{G(1/2)} = log{sqrt(Pi)} */

        if(lnbeta == 0) {
            lnbeta = LnGamma(df/2) + lnghalf - LnGamma((df+1)/2);
        }
        cdf = CDFBeta(df/(df+z*z), df/2, 0.5, lnbeta);

        if(z>=0) cdf = 1 - cdf/2;
        else     cdf /= 2;
        return(cdf);
    }

    private double CDFBeta (double x, double pin, double qin, double lnbeta)
    {
/* Returns distribution function of the standard form of the beta distribution,
   that is, the incomplete beta ratio I_x(p,q).

   This is also known as the incomplete beta function ratio I_x(p, q)

   lnbeta is log of the complete beta function; provide it if known,
   and otherwise use 0.

   This is called from QuantileBeta() in a root-finding loop.

    This routine is a translation into C of a Fortran subroutine
    by W. Fullerton of Los Alamos Scientific Laboratory.
    Bosten and Battiste (1974).
    Remark on Algorithm 179, CACM 17, p153, (1974).
*/
        double ans, c, finsum, p, ps, p1, q, term, xb, xi, y, small=1e-15;
        int n, i, ib;
        double eps = 0, alneps = 0, sml = 0, alnsml = 0;

        if(x<small)        return 0;
        else if(x>1-small) return 1;
        if(pin<=0 || qin<=0)  {
            throw new RuntimeException("pin qin parameter outside range in CDFBeta");
        }

        final double FLT_RADIX = 2;
        final double DBL_MANT_DIG = 53;

        if (eps == 0) {/* initialize machine constants ONCE */
            eps = Math.pow((double)FLT_RADIX, -(double)DBL_MANT_DIG);
            alneps = Math.log(eps);
            sml = Double.MIN_VALUE;
            alnsml = Math.log(sml);
        }
        y = x;  p = pin;  q = qin;

        /* swap tails if x is greater than the mean */
        if (p / (p + q) < x) {
            y = 1 - y;
            p = qin;
            q = pin;
        }

        if(lnbeta==0) lnbeta = LnBeta(p, q);

        if ((p + q) * y / (p + 1) < eps) {  /* tail approximation */
            ans = 0;
            xb = p * Math.log(Math.max(y, sml)) - Math.log(p) - lnbeta;
            if (xb > alnsml && y != 0)
                ans = Math.exp(xb);
            if (y != x || p != pin)
                ans = 1 - ans;
        }
        else {
            /* evaluate the infinite sum first.  term will equal */
            /* y^p / beta(ps, p) * (1 - ps)-sub-i * y^i / fac(i) */
            ps = q - Math.floor(q);
            if (ps == 0)
                ps = 1;

            xb=LnGamma(ps)+LnGamma(p)-LnGamma(ps+p);
            xb = p * Math.log(y) - xb - Math.log(p);

            ans = 0;
            if (xb >= alnsml) {
                ans = Math.exp(xb);
                term = ans * p;
                if (ps != 1) {
                    n = (int)Math.max(alneps/Math.log(y), 4.0);
                    for(i=1 ; i<= n ; i++) {
                        xi = i;
                        term = term * (xi - ps) * y / xi;
                        ans = ans + term / (p + xi);
                    }
                }
            }

            /* evaluate the finite sum. */
            if (q > 1) {
                xb = p * Math.log(y) + q * Math.log(1 - y) - lnbeta - Math.log(q);
                ib = (int) (xb/alnsml);  if(ib<0) ib=0;
                term = Math.exp(xb - ib * alnsml);
                c = 1 / (1 - y);
                p1 = q * c / (p + q - 1);

                finsum = 0;
                n = (int) q;
                if (q == (double)n)
                    n = n - 1;
                for(i=1 ; i<=n ; i++) {
                    if (p1 <= 1 && term / eps <= finsum)
                        break;
                    xi = i;
                    term = (q - xi + 1) * c * term / (p + q - xi);
                    if (term > 1) {
                        ib = ib - 1;
                        term = term * sml;
                    }
                    if (ib == 0)
                        finsum = finsum + term;
                }
                ans = ans + finsum;
            }
            if (y != x || p != pin)
                ans = 1 - ans;
            if(ans>1) ans=1;
            if(ans<0) ans=0;
        }
        return ans;
    }

    private double LnBeta(double p, double q) {
        return LnGamma(p) + LnGamma(q) - LnGamma(p + q);
    }

    @Override
    public double logPdf(double x) {
        return Math.log(pdf(x));
    }

    @Override
    public double cdf(double x) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double quantile(double y) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double mean() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double variance() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public UnivariateFunction getProbabilityDensityFunction() {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getInitialDomain(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getDomainLowerBound(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    protected double getDomainUpperBound(double v) {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public double cumulativeProbability(double v) throws MathException {
        throw new RuntimeException("Not yet implemented!");
    }

    @Override
    public Citation.Category getCategory() {
        return Citation.Category.PRIOR_MODELS;
    }

    @Override
    public String getDescription() {
        return "Skew-T distribution.";
    }

    @Override
    public List<Citation> getCitations() {
        return Collections.singletonList(CITATION);
    }

    public static Citation CITATION = new Citation(
            new Author[]{
                    new Author("RD", "Wilkingson"),
                    new Author("ME", "Steiper"),
                    new Author("C", "Soligo"),
                    new Author("RD", "Martin"),
                    new Author("Z", "Yang"),
                    new Author("S", "Tavaré")
            },
            "Dating Primate Divergences through an Integrated Analysis of Palaeontological and Molecular Data",
            2011,
            "Systematic Biology",
            60,
            16,
            31,
            Citation.Status.PUBLISHED
    );
}
