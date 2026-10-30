/*
 * CompactGradientTest.java
 *
 * Copyright (c) 2002-2026 Alexei Drummond, Andrew Rambaut and Marc Suchard
 */

package test.dr.inference.hmc;

import dr.inference.hmc.CompactGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.CompoundParameter;
import dr.inference.model.DuplicatedParameter;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import junit.framework.TestCase;

/** By the chain rule, compacting a gradient sums the entries of every dimension that maps onto one free value. */
public class CompactGradientTest extends TestCase {

    public void testDuplicatedParameterSumsCopies() {
        Parameter underlying = new Parameter.Default("underlying", new double[]{1.5, 2.0});
        DuplicatedParameter duplicated = duplicate(underlying, 3);

        CompactGradient compact = new CompactGradient(provider(duplicated));

        assertEquals(2, compact.getDimension());
        assertEquals(1.5, compact.getParameter().getParameterValue(0), 0.0);
        assertEquals(2.0, compact.getParameter().getParameterValue(1), 0.0);
        // copy-major layout: entries 0, 2, 4 belong to underlying[0] and entries 1, 3, 5 to underlying[1]
        assertArrayEquals(new double[]{1.0 + 3.0 + 5.0, 2.0 + 4.0 + 6.0}, compact.getGradientLogDensity());

        compact.getParameter().setParameterValue(1, 2.5);
        assertEquals(2.5, underlying.getParameterValue(1), 0.0);
        assertEquals(2.5, duplicated.getParameterValue(5), 0.0);
    }

    public void testDuplicatedParameterInsideCompoundParameter() {
        Parameter single = new Parameter.Default("single", 0.5);
        Parameter underlying = new Parameter.Default("underlying", 1.5);
        CompoundParameter compound = new CompoundParameter("compound",
                new Parameter[]{single, duplicate(underlying, 3), underlying});

        CompactGradient compact = new CompactGradient(provider(compound));

        assertEquals(2, compact.getDimension());
        assertArrayEquals(new double[]{1.0, 2.0 + 3.0 + 4.0 + 5.0}, compact.getGradientLogDensity());
    }

    public void testDuplicatedParameterOfCompoundAndDuplicatedParameters() {
        Parameter a = bounded(new Parameter.Default("a", new double[]{1.5, 2.0}));
        Parameter b = bounded(new Parameter.Default("b", 0.5));
        DuplicatedParameter ofCompound = duplicate(new CompoundParameter("ab", new Parameter[]{a, b}), 2);

        CompactGradient compact = new CompactGradient(provider(ofCompound));

        assertEquals(3, compact.getDimension());
        assertArrayEquals(new double[]{1.0 + 4.0, 2.0 + 5.0, 3.0 + 6.0}, compact.getGradientLogDensity());

        Parameter c = new Parameter.Default("c", new double[]{1.5, 2.0});
        DuplicatedParameter ofDuplicated = duplicate(duplicate(c, 2), 3);

        compact = new CompactGradient(provider(ofDuplicated));

        assertEquals(2, compact.getDimension());
        assertArrayEquals(new double[]{1.0 + 3.0 + 5.0 + 7.0 + 9.0 + 11.0, 2.0 + 4.0 + 6.0 + 8.0 + 10.0 + 12.0},
                compact.getGradientLogDensity());
    }

    /** A parameter first seen after a repeated one maps to its offset in the compacted parameter. */
    public void testUniqueParameterAfterRepeats() {
        Parameter underlying = new Parameter.Default("underlying", new double[]{1.5, 2.0});
        Parameter single = new Parameter.Default("single", 0.5);
        Parameter last = new Parameter.Default("last", 0.25);
        CompoundParameter compound = new CompoundParameter("compound",
                new Parameter[]{duplicate(underlying, 2), single, underlying, last, single});

        CompactGradient compact = new CompactGradient(provider(compound));

        assertEquals(4, compact.getDimension());
        assertArrayEquals(new double[]{1.0 + 3.0 + 6.0, 2.0 + 4.0 + 7.0, 5.0 + 9.0, 8.0},
                compact.getGradientLogDensity());
    }

    private static DuplicatedParameter duplicate(Parameter parameter, int copies) {
        if (parameter instanceof Parameter.Default) {
            bounded(parameter);
        }
        DuplicatedParameter duplicated = new DuplicatedParameter(parameter);
        duplicated.addDuplicationParameter(new Parameter.Default((double) copies));
        return duplicated;
    }

    /** A DuplicatedParameter copies the bounds of the parameter it duplicates, so that parameter needs some. */
    private static Parameter bounded(Parameter parameter) {
        parameter.addBounds(new Parameter.DefaultBounds(Double.POSITIVE_INFINITY, Double.NEGATIVE_INFINITY,
                parameter.getDimension()));
        return parameter;
    }

    /** Returns the gradient 1, 2, ..., dimension, so each compacted entry identifies the inputs it summed. */
    private static GradientWrtParameterProvider provider(final Parameter parameter) {
        return new GradientWrtParameterProvider() {
            public Likelihood getLikelihood() {
                return null;
            }

            public Parameter getParameter() {
                return parameter;
            }

            public int getDimension() {
                return parameter.getDimension();
            }

            public double[] getGradientLogDensity() {
                double[] gradient = new double[parameter.getDimension()];
                for (int i = 0; i < gradient.length; i++) {
                    gradient[i] = i + 1;
                }
                return gradient;
            }
        };
    }

    private static void assertArrayEquals(double[] expected, double[] actual) {
        assertEquals(expected.length, actual.length);
        for (int i = 0; i < expected.length; i++) {
            assertEquals("entry " + i, expected[i], actual[i], 0.0);
        }
    }
}
