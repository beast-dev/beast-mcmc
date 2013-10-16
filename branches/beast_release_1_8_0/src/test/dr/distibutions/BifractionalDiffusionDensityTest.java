package test.dr.distibutions;

import junit.framework.TestCase;
import dr.math.distributions.BifractionalDiffusionDensity;
import dr.math.distributions.NormalDistribution;

/**
 * @author Marc A. Suchard
 */
public class BifractionalDiffusionDensityTest extends TestCase {


       public void setUp() {

        // Do nothing
    }

    public void testBrownianPdf() {

        System.out.println("Testing normal density");
        double alpha = 2.0;
        double beta = 1.0;
        BifractionalDiffusionDensity density = new BifractionalDiffusionDensity(alpha, beta);

        double vLower = 1.0;
        double vUpper = 100.0;
        double vIncr = 10.0;

        double xLower = -5;
        double xUpper = +5;
        double xIncr = 1;

        for (double v = vLower; v <= vUpper; v *= vIncr) {
            for (double x = xLower; x <= xUpper; x += xIncr) {
                double pdf = density.pdf(x, v);
                double pdfCheck = NormalDistribution.pdf(x, 0, Math.sqrt(v));
                System.err.println("A = "+pdf+ "(x = "+x+", v = "+v+")");
                System.err.println("B = "+pdfCheck);
                assertEquals(pdf, pdfCheck, 1E-10);
            }
        }

        

    }



}
