package test.dr.inference.hmc;

import dr.inference.hmc.CompoundGradient;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.hmc.ScoreCompoundGradient;
import dr.inference.hmc.ScoreContributionsWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.inference.model.Parameter;
import dr.inferencexml.hmc.ScoreCompoundGradientParser;
import junit.framework.TestCase;

import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;

public class ScoreCompoundGradientTest extends TestCase {

    private static Likelihood likelihood() {
        return new Likelihood.Abstract(null) {
            @Override
            protected boolean getLikelihoodKnown() {
                return false;
            }

            @Override
            protected double calculateLogLikelihood() {
                return 0.0;
            }
        };
    }

    private static GradientWrtParameterProvider plain(final String name, final int dimension) {
        final Parameter parameter = new Parameter.Default(name, new double[dimension]);
        final Likelihood likelihood = likelihood();
        return new GradientWrtParameterProvider() {
            public Likelihood getLikelihood() { return likelihood; }
            public Parameter getParameter() { return parameter; }
            public int getDimension() { return dimension; }
            public double[] getGradientLogDensity() { return new double[dimension]; }
        };
    }

    private static GradientWrtParameterProvider scored(final String name, final double[][] rows) {
        final int dimension = rows[0].length;
        final Parameter parameter = new Parameter.Default(name, new double[dimension]);
        final Likelihood likelihood = likelihood();
        return new ScoreContributionsWrtParameterProvider() {
            public Likelihood getLikelihood() { return likelihood; }
            public Parameter getParameter() { return parameter; }
            public int getDimension() { return dimension; }
            public double[] getGradientLogDensity() { return new double[dimension]; }
            public double[][] getScoreContributionsLogDensity() { return rows; }
        };
    }

    public void testConcatenatesScoreRowsAndZeroFillsComponentsWithoutScores() {
        final List<GradientWrtParameterProvider> parts = new ArrayList<GradientWrtParameterProvider>();
        parts.add(scored("a", new double[][]{{1.0, 2.0}, {3.0, 4.0}}));
        parts.add(plain("b", 1));
        parts.add(scored("c", new double[][]{{5.0}, {6.0}}));
        final double[][] rows = new ScoreCompoundGradient(parts).getScoreContributionsLogDensity();
        assertEquals(2, rows.length);
        assertEquals(4, rows[0].length);
        assertEquals(1.0, rows[0][0], 0.0);
        assertEquals(2.0, rows[0][1], 0.0);
        assertEquals(0.0, rows[0][2], 0.0);
        assertEquals(5.0, rows[0][3], 0.0);
        assertEquals(3.0, rows[1][0], 0.0);
        assertEquals(6.0, rows[1][3], 0.0);
    }

    public void testNoScoreProvidersGivesNull() {
        final List<GradientWrtParameterProvider> parts = new ArrayList<GradientWrtParameterProvider>();
        parts.add(plain("a", 2));
        parts.add(plain("b", 1));
        assertNull(new ScoreCompoundGradient(parts).getScoreContributionsLogDensity());
    }

    public void testUnequalRowCountsAreRejected() {
        final List<GradientWrtParameterProvider> parts = new ArrayList<GradientWrtParameterProvider>();
        parts.add(scored("a", new double[][]{{1.0}, {2.0}}));
        parts.add(scored("b", new double[][]{{1.0}}));
        try {
            new ScoreCompoundGradient(parts).getScoreContributionsLogDensity();
            fail("Expected IllegalStateException");
        } catch (IllegalStateException expected) {
            // expected
        }
    }

    public void testPlainCompoundGradientDoesNotExposeScores() {
        assertFalse(ScoreContributionsWrtParameterProvider.class.isAssignableFrom(CompoundGradient.class));
    }

    public void testParserIsRegisteredForDevelopmentUse() throws Exception {
        assertEquals("scoreCompoundGradient", new ScoreCompoundGradientParser().getParserName());
        assertEquals(ScoreCompoundGradient.class, new ScoreCompoundGradientParser().getReturnType());
        final InputStream in = ScoreCompoundGradient.class.getResourceAsStream(
                "/dr/app/beast/development_parsers.properties");
        assertNotNull(in);
        final StringBuilder text = new StringBuilder();
        final byte[] buffer = new byte[4096];
        int n;
        while ((n = in.read(buffer)) > 0) {
            text.append(new String(buffer, 0, n, "ISO-8859-1"));
        }
        in.close();
        assertTrue(text.indexOf("dr.inferencexml.hmc.ScoreCompoundGradientParser") >= 0);
    }
}
