package dr.math;


import java.util.ArrayList;
import java.util.List;


class BSpline {


    public static double[] polyAdd(double[] p, double[] q) {
        int n = Math.max(p.length, q.length);
        double[] r = new double[n];
        for (int i = 0; i < n; i++) {
            double a = (i < p.length) ? p[i] : 0.0;
            double b = (i < q.length) ? q[i] : 0.0;
            r[i] = a + b;
        }
        return r;
    }

    public static double[] polyScale(double[] p, double s) {
        double[] r = new double[p.length];
        for (int i = 0; i < p.length; i++) {
            r[i] = p[i] * s;
        }
        return r;
    }

    public static double[] polyMultiply(double[] p, double[] q) {
        double[] r = new double[p.length + q.length - 1];
        for (int i = 0; i < p.length; i++) {
            for (int j = 0; j < q.length; j++) {
                r[i + j] += p[i] * q[j];
            }
        }
        return r;
    }

    public static double[] polySquare(double[] p) {
        return polyMultiply(p, p);
    }

    public static double polyIntegral(double[] p, double L, double U) {
        double total = 0.0;
        double uPow = U, lPow = L;
        for (int k = 0; k < p.length; k++) {
            total += p[k] / (k + 1) * (uPow - lPow);
            if (k < p.length - 1) {
                uPow *= U;
                lPow *= L;
            }
        }
        return total;
    }


    public static class PPoly {
        public final double[] knots;
        public final double[][] pieces;

        public PPoly(double[] knots, double[][] pieces) {
            this.knots = knots.clone();
            this.pieces = pieces;
        }

        public static PPoly zero(double[] knots) {
            int nIntervals = knots.length - 1;
            double[][] pieces = new double[nIntervals][];
            for (int i = 0; i < nIntervals; i++) {
                pieces[i] = new double[]{0.0};
            }
            return new PPoly(knots, pieces);
        }

        public PPoly scale(double s) {
            double[][] newPieces = new double[pieces.length][];
            for (int i = 0; i < pieces.length; i++) {
                newPieces[i] = polyScale(pieces[i], s);
            }
            return new PPoly(knots, newPieces);
        }

        public static PPoly add(PPoly a, PPoly b) {
            assert java.util.Arrays.equals(a.knots, b.knots);
            double[][] newPieces = new double[a.pieces.length][];
            for (int i = 0; i < a.pieces.length; i++) {
                newPieces[i] = polyAdd(a.pieces[i], b.pieces[i]);
            }
            return new PPoly(a.knots, newPieces);
        }

        public PPoly square() {
            double[][] newPieces = new double[pieces.length][];
            for (int i = 0; i < pieces.length; i++) {
                newPieces[i] = polySquare(pieces[i]);
            }
            return new PPoly(knots, newPieces);
        }

        public double integral(double a, double b) {
            double total = 0.0;
            for (int i = 0; i < knots.length - 1; i++) {
                double left = knots[i];
                double right = knots[i + 1];
                if (right <= a || left >= b) continue;
                double L = Math.max(a, left);
                double R = Math.min(b, right);
                if (L >= R) continue;
                total += polyIntegral(pieces[i], L, R);
            }
            return total;
        }


        public double evaluate(double x) {
            for (int i = 0; i < knots.length - 1; i++) {
                if (x >= knots[i] && (x < knots[i + 1] || i == knots.length - 2)) {
                    double[] p = pieces[i];
                    double val = 0.0;
                    double xPow = 1.0;
                    for (int k = 0; k < p.length; k++) {
                        val += p[k] * xPow;
                        xPow *= x;
                    }
                    return val;
                }
            }
            return 0.0;
        }

    }


    public static List<PPoly> bSplineBasis(double[] knots, int p) {
        int m = knots.length;
        int nIntervals = m - 1;

        List<PPoly> B_prev = new ArrayList<>();
        for (int i = 0; i < m - 1; i++) {
            double[][] pieces = new double[nIntervals][];
            for (int s = 0; s < nIntervals; s++) {
                pieces[s] = new double[]{0.0};
            }
            if (knots[i] < knots[i + 1]) {
                pieces[i] = new double[]{1.0};
            }
            B_prev.add(new PPoly(knots, pieces));
        }


        for (int k = 1; k <= p; k++) {
            List<PPoly> B_curr = new ArrayList<>();
            int nk = m - k - 1;
            for (int i = 0; i < nk; i++) {
                double denom1 = knots[i + k] - knots[i];
                double denom2 = knots[i + k + 1] - knots[i + 1];

                double[][] newPieces = new double[nIntervals][];
                for (int s = 0; s < nIntervals; s++) {
                    double[] poly = {0.0};

                    if (denom1 != 0.0) {
                        double[] lin1 = {-knots[i], 1.0};
                        double[] prod = polyMultiply(lin1, B_prev.get(i).pieces[s]);
                        double[] scaled = polyScale(prod, 1.0 / denom1);
                        poly = polyAdd(poly, scaled);
                    }
                    if (denom2 != 0.0) {
                        double[] lin2 = {knots[i + k + 1], -1.0};
                        double[] prod = polyMultiply(lin2, B_prev.get(i + 1).pieces[s]);
                        double[] scaled = polyScale(prod, 1.0 / denom2);
                        poly = polyAdd(poly, scaled);
                    }
                    newPieces[s] = poly;
                }
                B_curr.add(new PPoly(knots, newPieces));
            }
            B_prev = B_curr;
        }
        return B_prev;
    }

}