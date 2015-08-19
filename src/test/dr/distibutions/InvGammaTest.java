package test.dr.distibutions;

import junit.framework.TestCase;
import dr.math.distributions.InverseGammaDistribution;

/**
 * Simple test for inverse gamma distribution.
 *
 * @author Joseph Heled
 *         Date: 24/04/2009
 */
public class InvGammaTest extends TestCase {
    interface TestData {
        double getShape();

        double getScale();

        double[] getPDF();

        double[] getCDF();
    }

    // test data generated from this python code:
// import scipy.stats
//
//print """TestData[] tests = {"""
//for shape,scale in ((3,2), (3,1)) :
//  d = scipy.stats.invgamma(shape, scale=scale)
//  print """
//        new TestData() {
//            public double getShape() {
//                return %d;
//            }
//
//            public double getScale() {
//               return %d;
//            }""" % (shape, scale)
//  x = (0.5, 1, 2)
//  print """
//            public double[] getPDF() {
//                return new double[]{%s};
//            }""" % " , ".join(["%g,%.14lf" % (z,d.pdf(z)) for z in x])
//
//  print """
//            public double[] getCDF() {
//                return new double[]{%s};
//            }
//        } ,""" % " , ".join(["%g,%.14lf" % (z,d.cdf(z)) for z in x])
//

    TestData[] tests = {

            new TestData() {
                public double getShape() {
                    return 3;
                }

                public double getScale() {
                    return 2;
                }

                public double[] getPDF() {
                    return new double[]{0.5, 1.17220088887899, 1, 0.54134113294645, 2, 0.09196986029286};
                }

                public double[] getCDF() {
                    return new double[]{0.5, 0.23810330555354, 1, 0.67667641618306, 2, 0.91969860292861};
                }
            },

            new TestData() {
                public double getShape() {
                    return 3;
                }

                public double getScale() {
                    return 1;
                }

                public double[] getPDF() {
                    return new double[]{0.5, 1.08268226589290, 1, 0.18393972058572, 2, 0.01895408311602};
                }

                public double[] getCDF() {
                    return new double[]{0.5, 0.67667641618306, 1, 0.91969860292861, 2, 0.98561232203303};
                }
            },
    };


    public void testInvGamma() {
        for( TestData td : tests ) {
            InverseGammaDistribution d = new InverseGammaDistribution(td.getShape(), td.getScale());

            {
                double[] p = td.getPDF();
                for(int k = 0; k < p.length; k += 2) {
                    assertEquals(d.pdf(p[k]), p[k + 1], 1e-10);

                    assertEquals(d.logPdf(p[k]), Math.log(p[k + 1]), 1e-10);
                }
            }

            double[] cdf = td.getCDF();
            for(int k = 0; k < cdf.length; k += 2) {
                assertEquals(d.cdf(cdf[k]), cdf[k + 1], 1e-10);
            }

//            int count[] = new int[cdf.length];
//            final int N = 100000;
//            for(int k = 0; k < N; ++k) {
//                double x = d.nextInverseGamma();
//                for(int l = 0; l < cdf.length; l += 2) {
//                    if( x < cdf[l] ) {
//                        count[l / 2] += 1;
//                    }
//                }
//            }
//            for(int l = 0; l < cdf.length; l += 2) {
//                assertEquals(count[l / 2] / (double) N, cdf[l + 1], 5e-3);
//            }
        }
    }
}
