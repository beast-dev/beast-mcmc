/*
 * Polynomial.java
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

package dr.math;

//import org.apfloat.Apfloat;
//import org.apfloat.ApfloatMath;

import java.math.BigDecimal;
import java.math.MathContext;
import java.util.Arrays;

/**
 * @author Marc A. Suchard
 */

public interface Polynomial extends Cloneable {

    public int getDegree();

    public Polynomial multiply(Polynomial b);

    public Polynomial integrate();

    public double evaluate(double x);

    public double logEvaluate(double x);

    public double logEvaluateHorner(double x);

    public void expand(double x);

    public Polynomial integrateWithLowerBound(double bound);

    public double getCoefficient(int n);

    public String getCoefficientString(int n);

    public void setCoefficient(int n, double x);

    public Polynomial getPolynomial();

    public Polynomial copy();

    public abstract class Abstract implements Polynomial {

        public abstract int getDegree();

        public abstract Polynomial multiply(Polynomial b);

        public abstract Polynomial integrate();

        public abstract double evaluate(double x);

        public abstract double getCoefficient(int n);

        public abstract void setCoefficient(int n, double x);

        public abstract Polynomial integrateWithLowerBound(double bound);

        public Polynomial getPolynomial() { return this; }

        public abstract double logEvaluate(double x);

        public abstract double logEvaluateHorner(double x);

        public abstract void expand(double x);

        public String toString() {
            StringBuffer bf = new StringBuffer();
            for (int n = getDegree(); n >= 0; n--) {
                bf.append(getCoefficientString(n));
                bf.append(X);
                bf.append(n);
                if (n > 0)
                    bf.append(" + ");
            }
            return bf.toString();
        }

        public abstract String getCoefficientString(int n);


        protected static final String FORMAT = "%3.2e";
        private static final String X = " x^";
    }

    public class LogDouble extends Abstract {

        public LogDouble(double[] coefficient) {
            this.logCoefficient = new double[coefficient.length];
            this.positiveCoefficient = new boolean[coefficient.length];
            for(int i=0; i<coefficient.length; i++) {
                if(coefficient[i] < 0) {
                    this.logCoefficient[i] = Math.log(-coefficient[i]);
                    this.positiveCoefficient[i] = false;
                } else {
                    this.logCoefficient[i] = Math.log(coefficient[i]);
                    this.positiveCoefficient[i] = true;
                }
            }
        }

        public double getLogCoefficient(int n) {
            return logCoefficient[n];
        }

        public void expand(double x) {
            final int degree = getDegree();
            for(int i=0; i<=degree; i++)
                logCoefficient[i] = x + logCoefficient[i];
        }

        public String getCoefficientString(int n) {
                 return String.format(FORMAT, getCoefficient(n));
             }


        public LogDouble(double[] logCoefficient, boolean[] positiveCoefficient) {
            this.logCoefficient = logCoefficient;
            if (positiveCoefficient != null)
                this.positiveCoefficient = positiveCoefficient;
            else {
                this.positiveCoefficient = new boolean[logCoefficient.length];
                Arrays.fill(this.positiveCoefficient,true);
            }
        }

        public LogDouble copy() {
            return new LogDouble(logCoefficient.clone(), positiveCoefficient.clone());
        }

        public int getDegree() {
            return logCoefficient.length - 1;
        }

        public LogDouble multiply(Polynomial inB) {
            if (!(inB.getPolynomial() instanceof LogDouble))
                throw new RuntimeException("yuck!");

            LogDouble b = (LogDouble) inB.getPolynomial();
            final int degreeA = getDegree();
            final int degreeB = b.getDegree();

            double[] newLogCoefficient = new double[degreeA + degreeB + 1];
            boolean[] newPositiveCoefficient = new boolean[degreeA + degreeB + 1];
            Arrays.fill(newLogCoefficient, java.lang.Double.NEGATIVE_INFINITY);
            Arrays.fill(newPositiveCoefficient, true);

            for (int n = 0; n <= degreeA; n++) {
                for (int m = 0; m <= degreeB; m++) {
                    final double change = logCoefficient[n] + b.logCoefficient[m];
                    final int nm = n + m;
                    final boolean positiveChange = !(positiveCoefficient[n] ^ b.positiveCoefficient[m]);
                    if (newLogCoefficient[nm] == java.lang.Double.NEGATIVE_INFINITY) {
                        newLogCoefficient[nm] = change;
                        newPositiveCoefficient[nm] = positiveChange;
                    } else {
                        if (change != 0.0) {
                            if (newPositiveCoefficient[nm] ^ positiveChange) { // Sign difference, must subtract
                                if (newLogCoefficient[nm] > change)
                                    newLogCoefficient[nm] = LogTricks.logDiff(newLogCoefficient[nm], change);
                                else {
                                    newLogCoefficient[nm] = LogTricks.logDiff(change, newLogCoefficient[nm]);
                                    newPositiveCoefficient[nm] = !newPositiveCoefficient[nm]; // Switch signs
                                }
                            } else { // Same signs, just add
                                newLogCoefficient[nm] = LogTricks.logSum(newLogCoefficient[nm], change);
                            }
                        }
                    }
                }
            }
            return new LogDouble(newLogCoefficient,newPositiveCoefficient);
        }

        public LogDouble integrate() {
            final int degree = getDegree();
            double[] newLogCoefficient = new double[degree + 2];
            boolean[] newPositiveCoefficient = new boolean[degree + 2];
            for (int n=0; n<=degree; n++) {
                newLogCoefficient[n+1] = logCoefficient[n] - Math.log(n+1);
                newPositiveCoefficient[n+1] = positiveCoefficient[n];
            }
            newLogCoefficient[0] = java.lang.Double.NEGATIVE_INFINITY;
            newPositiveCoefficient[0] = true;
            return new LogDouble(newLogCoefficient,newPositiveCoefficient);
        }

        public double evaluate(double x) {
            SignedLogDouble result = signedLogEvaluate(x);
            double value = Math.exp(result.value);
            if (!result.positive)
                value = -value;
            return value;
        }

        public double evaluateAsReal(double x) {
            double result = 0;
            double xn = 1;
            for (int n = 0; n <= getDegree(); n++) {
                result += xn * getCoefficient(n);
                xn *= x;
            }
            return result;

        }

        public double logEvaluate(double x) {
            if (x < 0)
                throw new RuntimeException("Negative arguments not yet implemented in Polynomial.LogDouble");
            SignedLogDouble result = signedLogEvaluate(x);
            if (result.positive)
                return result.value;
            return java.lang.Double.NaN;
//            return -result.value;
        }

        public double logEvaluateHorner(double x) {
            if (x < 0)
                throw new RuntimeException("Negative arguments not yet implemented in Polynomial.LogDouble");
            SignedLogDouble result = signedLogEvaluateHorners(x);
            if (result.positive)
                return result.value;
            return java.lang.Double.NaN;
//            return -result.value;
        }


        public SignedLogDouble signedLogEvaluateHorners(double x) {
            // Using Horner's Rule
            final double logX = Math.log(x);
            final int degree = getDegree();
            double logResult = logCoefficient[degree];
            boolean positive = positiveCoefficient[degree];
            for(int n=degree-1; n>=0; n--) {
                logResult += logX;
                if (!(positiveCoefficient[n] ^ positive)) // Same sign
                    logResult = LogTricks.logSum(logResult,logCoefficient[n]);
                else { // Different signs
                    if (logResult > logCoefficient[n])
                        logResult = LogTricks.logDiff(logResult,logCoefficient[n]);
                    else {
                        logResult = LogTricks.logDiff(logCoefficient[n],logResult);
                        positive = !positive;
                    }
                }
            }
            return new SignedLogDouble(logResult,positive);
        }

        private SignedLogDouble signedLogEvaluate(double x) {
            final double logX = Math.log(x);
            final int degree = getDegree();
            double logResult = logCoefficient[0];
            boolean positive = positiveCoefficient[0];
            for(int n=1; n<=degree; n++) {
//                logResult += logX;
                final double value = n*logX + logCoefficient[n];
                if (!(positiveCoefficient[n] ^ positive)) // Same sign
                    logResult = LogTricks.logSum(logResult,value);
                else { // Different signs
                    if (logResult > value)
                        logResult = LogTricks.logDiff(logResult,value);
                    else {
                        logResult = LogTricks.logDiff(value,logResult);
                        positive = !positive;
                    }
                }
            }
            return new SignedLogDouble(logResult,positive);
        }

        public double getCoefficient(int n) {
            double coef = Math.exp(logCoefficient[n]);
            if (!positiveCoefficient[n])
                coef *= -1;
            return coef;
        }

        public void setCoefficient(int n, double x) {
            if (x < 0) {
                positiveCoefficient[n] = false;
                x = -x;
            } else
                positiveCoefficient[n] = true;
            logCoefficient[n] = Math.log(x);
        }

        public Polynomial integrateWithLowerBound(double bound) {
            LogDouble integrand = integrate();
            SignedLogDouble signedLogDouble = integrand.signedLogEvaluate(bound);
            integrand.logCoefficient[0] = signedLogDouble.value;
            integrand.positiveCoefficient[0] = !signedLogDouble.positive;
            return integrand;
        }

        double[] logCoefficient;
        boolean[] positiveCoefficient;

        class SignedLogDouble {

            double value;
            boolean positive;

            SignedLogDouble(double value, boolean positive) {
                this.value = value;
                this.positive = positive;
            }
        }
    }

    public class BigDouble extends Abstract {

        private static MathContext precision = new MathContext(1000);

        public BigDouble(double[] doubleCoefficient) {
            this.coefficient = new BigDecimal[doubleCoefficient.length];
            for(int i=0; i<doubleCoefficient.length; i++)
                coefficient[i] = new BigDecimal(doubleCoefficient[i]);
        }

        public BigDouble copy() {
            return new BigDouble(coefficient.clone());
        }

        public String getCoefficientString(int n) {
            return coefficient[n].toString();

             }

         public void expand(double x) {
            throw new RuntimeException("Not yet implement: Polynomial.BigDouble.expand()");
        }



        public BigDouble(BigDecimal[] coefficient) {
            this.coefficient = coefficient;
        }

        public int getDegree() {
             return coefficient.length - 1;
        }

        public BigDouble multiply(Polynomial b) {
            if (!(b.getPolynomial() instanceof BigDouble))
                throw new RuntimeException("Incompatiable polynomial types");
            BigDouble bd = (BigDouble) b.getPolynomial();
            BigDecimal[] newCoefficient = new BigDecimal[getDegree() + bd.getDegree()+1];
            for(int i=0; i<newCoefficient.length; i++)
                newCoefficient[i] = new BigDecimal(0.0);
            for(int n=0; n<=getDegree(); n++) {
                for(int m=0; m<=bd.getDegree(); m++)
                   newCoefficient[n+m] = newCoefficient[n+m].add(coefficient[n].multiply(bd.coefficient[m]));
            }
            return new BigDouble(newCoefficient);
        }

        public BigDouble integrate() {
            BigDecimal[] newCoefficient = new BigDecimal[getDegree()+2];
            for(int n=0; n<=getDegree(); n++) {
                newCoefficient[n+1] = coefficient[n].divide(new BigDecimal(n+1),precision);
            }
            newCoefficient[0] = new BigDecimal(0.0);
            return new BigDouble(newCoefficient);
        }

        public double evaluate(double x) {
            return evaluateBigDecimal(new BigDecimal(x)).doubleValue();
        }

        public double logEvaluate(double x) {
            return BigDecimalUtils.ln(evaluateBigDecimal(new BigDecimal(x)),10).doubleValue();
        }

        public double logEvaluateHorner(double x) {
            return logEvaluate(x);
        }

        protected BigDecimal evaluateBigDecimal(BigDecimal x) {
            BigDecimal result = new BigDecimal(0.0);
            BigDecimal xn = new BigDecimal(1.0);
            for(int n=0; n<=getDegree(); n++) {
                result = result.add(coefficient[n].multiply(xn));
                xn = xn.multiply(x);
            }
            return result;
        }

        public double getCoefficient(int n) {
            return coefficient[n].doubleValue();
        }

        public void setCoefficient(int n, double x) {
            coefficient[n] = new BigDecimal(x);
        }

        public Polynomial integrateWithLowerBound(double bound) {
            BigDouble integrand = integrate();
            final BigDecimal x0 = integrand.evaluateBigDecimal(new BigDecimal(bound));
            integrand.coefficient[0] = x0.multiply(new BigDecimal(-1.0));
            return integrand;
        }

        public void setCoefficient(int n, BigDecimal x) {
            coefficient[n] = x;
        }

        BigDecimal[] coefficient;
    }

//    public class APDouble extends Abstract {
//
//        public String getCoefficientString(int n) {
//            return coefficient[n].toString();
//        }
//
//        public APDouble copy() {
//            return new APDouble(coefficient.clone());
//        }
//
//         public APDouble(double[] doubleCoefficient) {
//             this.coefficient = new Apfloat[doubleCoefficient.length];
//             for(int i=0; i<doubleCoefficient.length; i++)
//                 coefficient[i] = new Apfloat(doubleCoefficient[i]);
//         }
//
//         public APDouble(Apfloat[] coefficient) {
//             this.coefficient = coefficient;
//         }
//
//         public int getDegree() {
//              return coefficient.length - 1;
//         }
//
//         public APDouble multiply(Polynomial b) {
//             if (!(b.getPolynomial() instanceof APDouble))
//                 throw new RuntimeException("Incompatiable polynomial types");
//             APDouble bd = (APDouble) b.getPolynomial();
//             Apfloat[] newCoefficient = new Apfloat[getDegree() + bd.getDegree()+1];
//             for(int i=0; i<newCoefficient.length; i++)
//                 newCoefficient[i] = new Apfloat(0.0);
//             for(int n=0; n<=getDegree(); n++) {
//                 for(int m=0; m<=bd.getDegree(); m++)
//                    newCoefficient[n+m] = newCoefficient[n+m].add(coefficient[n].multiply(bd.coefficient[m]));
//             }
//             return new APDouble(newCoefficient);
//         }
//
//         public APDouble integrate() {
//             Apfloat[] newCoefficient = new Apfloat[getDegree()+2];
//             for(int n=0; n<=getDegree(); n++) {
//                 newCoefficient[n+1] = coefficient[n].divide(new Apfloat(n+1));
//             }
//             newCoefficient[0] = new Apfloat(0.0);
//             return new APDouble(newCoefficient);
//         }
//
//         public double evaluate(double x) {
//             return evaluateAPDouble(new Apfloat(x)).doubleValue();
//         }
//
//         public double logEvaluate(double x) {
//             Apfloat result = evaluateAPDouble(new Apfloat((x)));
//             if (result.doubleValue() == 0)
//                return java.lang.Double.NEGATIVE_INFINITY;
//             return ApfloatMath.log(result).doubleValue();
//         }
//
//         public double logEvaluateHorner(double x) {
//             return logEvaluateInLogSpace(x);
//         }
//
//         private static double log(Apfloat x) {
//             double log = ApfloatMath.log(x).doubleValue();
//             if (java.lang.Double.isInfinite(log))
//                 throw new RuntimeException("Still infinite");
//             return log;
//         }
//
//        private static boolean positive(Apfloat x) {
//            return x.signum() != -1;
//        }
//
//         public double logEvaluateInLogSpace(double x) {
//             // Using Horner's Rule
//             final double logX = Math.log(x);
//             final int degree = getDegree();
//             boolean positive = positive(coefficient[degree]);
//             double logResult;
//             if (positive)
//                logResult = log(coefficient[degree]);
//             else
//                logResult = log(coefficient[degree].negate());
//             for(int n=degree-1; n>=0; n--) {
//                 logResult += logX;
//                 if (coefficient[n].signum() != 0) {
//                 final boolean nextPositive = positive(coefficient[n]);
//                 double logNextValue;
//                 if (nextPositive)
//                    logNextValue = log(coefficient[n]);
//                 else
//                    logNextValue = log(coefficient[n].negate());
//                 if(!(nextPositive ^ positive)) // Same sign
//                    logResult = LogTricks.logSum(logResult,logNextValue);
//                 else { // Different signs
//                     if (logResult > logNextValue)
//                         logResult = LogTricks.logDiff(logResult,logNextValue);
//                     else {
//                         logResult = LogTricks.logDiff(logNextValue,logResult);
//                         positive = !positive;
//                     }
//                 }
//                 }
//             }
//             if (!positive)
//                logResult = -logResult;
//             return logResult;
//         }
//
//         protected Apfloat evaluateAPDouble(Apfloat x) {
//             Apfloat result = new Apfloat(0.0);
//             Apfloat xn = new Apfloat(1.0);
//             for(int n=0; n<=getDegree(); n++) {
//                 result = result.add(coefficient[n].multiply(xn));
//                 xn = xn.multiply(x);
//             }
//             // TODO Rewrite using Horner's Rule
//             return result;
//         }
//
//         public double getCoefficient(int n) {
//             return coefficient[n].doubleValue();
//         }
//
//         public void setCoefficient(int n, double x) {
//             coefficient[n] = new Apfloat(x);
//         }
//
//         public Polynomial integrateWithLowerBound(double bound) {
//             APDouble integrand = integrate();
//             final Apfloat x0 = integrand.evaluateAPDouble(new Apfloat(bound));
//             integrand.coefficient[0] = x0.multiply(new Apfloat(-1.0));
//             return integrand;
//         }
//
//         public void setCoefficient(int n, Apfloat x) {
//             coefficient[n] = x;
//         }
//
//         Apfloat[] coefficient;
//     }


    public class Double extends Abstract {

        public Double(double[] coefficient) {
            this.coefficient = coefficient;
        }

        public Double copy() {
            return new Double(coefficient.clone());
        }

        public Double(Polynomial polynomial) {
            this.coefficient = new double[polynomial.getDegree() + 1];
            for (int n = 0; n <= polynomial.getDegree(); n++)
                coefficient[n] = polynomial.getCoefficient(n);
        }

        public void expand(double x) {
            final int degree = getDegree();
            for(int i=0; i<=degree; i++)
                coefficient[i] = x * coefficient[i];
        }

        public String getCoefficientString(int n) {
                 return String.format(FORMAT, getCoefficient(n));
             }


        public int getDegree() {
            return coefficient.length - 1;
        }

        public double getCoefficient(int n) {
            return coefficient[n];
        }

        public double logEvaluate(double x) {
            return Math.log(evaluate(x));
        }

        public double logEvaluateQuick(double x, int n) {
            return Math.log(evaluateQuick(x,n));
        }

        public double logEvaluateHorner(double x) {
            // Uses Horner's Rule in log scale
            final int degree = getDegree();
            final double logX = Math.log(x);
            boolean positive = coefficient[degree] > 0;
            double logResult;
            if (positive)
                logResult = Math.log(coefficient[degree]);
            else
                logResult = Math.log(-coefficient[degree]);
            for(int n=degree-1; n>=0; n--) {
                logResult += logX;
                boolean positiveCoefficient = coefficient[n] > 0;
                double logCoefficient;
                if (positiveCoefficient)
                    logCoefficient = Math.log(coefficient[n]);
                else
                    logCoefficient = Math.log(-coefficient[n]);
                if (!(positiveCoefficient ^ positive)) // Same sign
                    logResult = LogTricks.logSum(logResult,logCoefficient);
                else { // Different signs
                    if (logResult > logCoefficient)
                        logResult = LogTricks.logDiff(logResult,logCoefficient);
                    else {
                        logResult = LogTricks.logDiff(logCoefficient,logResult);
                        positive = !positive;
                    }
                }
            }
            if (!positive)
                return java.lang.Double.NaN;
            return logResult;
        }

        public Polynomial.Double multiply(Polynomial b) {
            double[] newCoefficient = new double[getDegree() + b.getDegree() + 1];
            for (int n = 0; n <= getDegree(); n++) {
                for (int m = 0; m <= b.getDegree(); m++) {
                    newCoefficient[n + m] += coefficient[n] * b.getCoefficient(m);
                }
            }
            return new Double(newCoefficient);
        }

        public Polynomial.Double integrate() {
            double[] newCoefficient = new double[getDegree() + 2];
            for (int n = 0; n <= getDegree(); n++) {
                newCoefficient[n + 1] = coefficient[n] / (n + 1);
            }
            return new Double(newCoefficient);
        }

        public double evaluateHorners(double x) {
            // Uses Horner's Rule
            final int degree = getDegree();
            double result = coefficient[degree];
            for (int n=degree-1; n>=0; n--)
                result = result*x + coefficient[n];
            return result;
        }

        public double evaluateQuick(double x, int n) {
            int m = getDegree();
            double xm = Math.pow(x,m);
            double result = xm * coefficient[m];
            for (int i=n-1; i>0; i--) {
                xm /= x;
                m--;
                result += xm * coefficient[m];
            }
            return result;
        }

        public double evaluate(double x) {
            double result = 0.0;
            double xn = 1;
            for (int n = 0; n <= getDegree(); n++) {
                result += xn * coefficient[n];
                xn *= x;
            }
            return result;
        }

        public Polynomial integrateWithLowerBound(double bound) {
            Double integrand = integrate();
//            System.err.println("integrand = "+integrand);
            integrand.coefficient[0] = -integrand.evaluate(bound);
            return integrand;
        }

        public void setCoefficient(int n, double x) {
            coefficient[n] = x;
        }

        double[] coefficient;
    }

    public enum Type {
        DOUBLE, APDOUBLE, LOG_DOUBLE, BIG_DOUBLE//, RATIONAL
    }
}
