package dr.inference.distribution;

import dr.math.MathUtils;
import dr.math.distributions.RandomGenerator;

public class ExponentialTiltedStableDistribution implements RandomGenerator {

    @Override
    public Object nextRandom() {
        return null;
    }

    @Override
    public double logPdf(Object x) {
        return 0;
    }

    private static final double doubleRejectionCost = 4.0;

    public static double nextTiltedStable(double exponent, double tilt) {
        double divideAndConquerCost = Math.pow(tilt, exponent);
        if (divideAndConquerCost < doubleRejectionCost) {
            return divideAndConquerSample(exponent, tilt);
        } else {
            return doubleRejectionSample(exponent, tilt);
        }
    }

    private static double divideAndConquerSample(double alpha, double lambda) {

        int partitionSize = Math.max(1, (int)Math.floor(Math.pow(lambda, alpha)));
        double c = Math.pow(1.0 / partitionSize, 1.0 / alpha);

        double x = 0.0;
        for (int i = 0; i < partitionSize; ++i) {
            x += dividedRvSample(alpha, lambda, c);
        }

        return x;
    }

    private static double dividedRvSample(double alpha, double lambda, double c) {
        boolean accept = false;
        double s = 0.0;
        while (!accept) {
            s = c * sampleNonTiltedRv(alpha);
            accept = MathUtils.nextDouble() < Math.exp(-lambda * s);
        }

        return s;
    }

    private static double sampleNonTiltedRv(double alpha) {
        double v = MathUtils.nextDouble();
        double e = -Math.log(MathUtils.nextDouble());
        return Math.pow(
                zolotarevFunction(Math.PI * v, alpha) / e,
                (1.0 - alpha) / alpha);
    }

    private static double doubleRejectionSample(double alpha, double lambda) {
        return retstable_LD(lambda, alpha);
    }

    private static double zolotarevFunction(double x, double alpha) {
        return Math.pow(
                Math.pow(
                        (1.0 - alpha) * sinc((1.0 - alpha) * x),
                        1.0 - alpha)
                        * Math.pow(
                        alpha * sinc(alpha * x),
                        alpha)
                        / sinc(x)
                , 1.0 / (1.0 - alpha)
        );
    }

    private static final double M_PI = Math.PI;
    private static final double M_SQRT_PI = Math.sqrt(Math.PI);
    private static final double M_SQRT2 = Math.sqrt(2.0);
    private static final double M_PI_2 = 0.5 * Math.PI;

    /*
      Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan

      This program is free software; you can redistribute it and/or modify it under
      the terms of the GNU General Public License as published by the Free Software
      Foundation; either version 3 of the License, or (at your option) any later
      version.
    */
    private static double retstable_LD(final double h, final double alpha) {
        
        /*
          alpha == 1 => St corresponds to a point mass at V0 with Laplace-Stieltjes
          transform exp(-V0*t)
         */
        if (alpha == 1.) {
            return 1.0;
        }

        // compute variables not depending on V0
        final double c1 = Math.sqrt(M_PI_2);
        final double c2 = 2. + c1;
        double Ialpha = 1. - alpha,
                b = Ialpha / alpha,
                lambda_alpha = Math.pow(h, alpha); // < Marius Hofert: work directly with lambda^alpha (numerically more stable for small alpha)

        /*
          Apply the algorithm of Devroye (2009) to draw from
          \tilde{S}(alpha, 1, (cos(alpha*pi/2))^{1/alpha}, I_{alpha = 1},
          	     lambda*I_{alpha != 1};1) with Laplace-Stieltjes transform
          exp(-((lambda+t)^alpha-lambda^alpha))
         */
        double gamma = lambda_alpha * alpha * Ialpha;
        double sgamma = Math.sqrt(gamma);
        double c3 = c2 * sgamma;
        double xi = (1. + M_SQRT2 * c3) / M_PI; /*< according to John Lau */
        double psi = c3 * Math.exp(-gamma * M_PI * M_PI / 8.) / M_SQRT_PI;
        double w1 = c1 * xi / sgamma;
        double w2 = 2. * M_SQRT_PI * psi;
        double w3 = xi * M_PI;

        double X, c, E;
        do {
            double U, z, Z;
            do {
                double V = MathUtils.nextDouble();
                if (gamma >= 1) {
                    if (V < w1 / (w1 + w2)) U = Math.abs(MathUtils.nextGaussian()) / sgamma;
                    else {
                        double W_ = MathUtils.nextDouble();
                        U = M_PI * (1. - W_ * W_);
                    }
                } else {
                    double W_ = MathUtils.nextDouble();
                    if (V < w3 / (w2 + w3)) U = M_PI * W_;
                    else U = M_PI * (1. - W_ * W_);
                }
                double W = MathUtils.nextDouble();
                double zeta = Math.sqrt(BdB0(U, alpha));
                z = 1 / (1 - Math.pow(1 + alpha * zeta / sgamma, -1 / alpha)); /*< Marius Hofert: numerically more stable for small alpha */
                /*< compute rho */
                double rho = M_PI * Math.exp(-lambda_alpha * (1. - 1. / (zeta * zeta))) /
                        ((1. + c1) * sgamma / zeta + z);
                double d = 0.;
                if (U >= 0 && gamma >= 1) d += xi * Math.exp(-gamma * U * U / 2.);
                if (U > 0 && U < M_PI) d += psi / Math.sqrt(M_PI - U);
                if (U >= 0 && U <= M_PI && gamma < 1) d += xi;
                rho *= d;
                Z = W * rho;
            } while (!(U < M_PI && Z <= 1.)); /* check rejection condition */

            double
                    a = zolotarevFunction(U, alpha), //Math.pow(_A_3(U,alpha,Ialpha), 1./Ialpha),
                    m = Math.pow(b / a, alpha) * lambda_alpha,
                    delta = Math.sqrt(m * alpha / a),
                    a1 = delta * c1,
                    a3 = z / a,
                    s = a1 + delta + a3;
            double V_ = MathUtils.nextDouble(), N_ = 0., E_ = 0.; // -Wall
            if (V_ < a1 / s) {
                N_ = MathUtils.nextGaussian();
                X = m - delta * Math.abs(N_);
            } else {
                if (V_ < (a1 + delta) / s) /*< according to John Lau */
                    X = m + delta * MathUtils.nextDouble();
                else {
                    E_ = -Math.log(MathUtils.nextDouble());
                    X = m + delta + E_ * a3;
                }
            }
            E = -Math.log(Z);
            /*< check rejection condition */
            c = a * (X - m) + Math.exp((1 / alpha) * Math.log(lambda_alpha) - b * Math.log(m)) * (Math.pow(m / X, b) - 1); /*< Marius Hofert: numerically more stable for small alpha */
            if (X < m) c -= N_ * N_ / 2.;
            else if (X > m + delta) c -= E_;

        } while (!(X >= 0 && c <= E));
        /*
          Transform variates from the distribution corresponding to the
          Laplace-Stieltjes transform exp(-((lambda+t)^alpha-lambda^alpha))
          to those of the distribution corresponding to the Laplace-Stieltjes
          transform exp(-V_0((h+t)^alpha-h^alpha)).
         */

        return Math.exp(-b * Math.log(X)); /*< Marius Hofert: numerically more stable for small alpha */
    }
    
    /**
     * Evaluation of B(x)/B(0), see Devroye (2009).
     *
     * @param x argument
     * @param alpha parameter in (0,1]
     * @return sinc(x) / (sinc(alpha*x)^alpha * sinc((1-alpha)*x)^(1-alpha))
     * @author Martin Maechler (2010-04-28)
     */
    private static double BdB0(double x, double alpha) {
        double Ialpha = 1. - alpha;
        double den = Math.pow(sinc(alpha * x), alpha) * Math.pow(sinc(Ialpha * x), Ialpha);
        return sinc(x) / den;
    }

    /**
     * Fast and accurate evaluation of sinc(x) := sin(x)/x, including the limit x=0
     *
     * @param x any (double precision) number
     * @return sinc(x)
     * @author Martin Maechler (2010-04-28)
     */
    private static double sinc(double x) {
        double ax = Math.abs(x);
        if (ax < 0.006) {
            if (x == 0.) { return 1; }
            double x2 = x * x;
            if (ax < 2e-4) {
                return 1. - x2 / 6.;
            } else {
                return 1. - x2 / 6. * (1 - x2 / 20.);
            }
        } else {
            return Math.sin(x) / x;
        }
    }

//    private static double sinc(double x) {
//        if (Math.abs(x) < 0.01) {
//            double xx = x * x;
//            return 1.0 - xx / 6.0 * (1.0 - xx / 20.0);
//        } else {
//            return Math.sin(x) / x;
//        }
}
