package dr.evomodelxml.continuous.hmc;

import dr.evomodel.treedatalikelihood.discrete.GradientBenchmark;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.xml.*;

/**
 * XML parser for {@link GradientBenchmark}.
 *
 * Syntax:
 * <pre>{@code
 * <benchmarkGradient warmup="50" iterations="200">
 *   <!-- any GradientWrtParameterProvider, e.g. exactLogCtmcRateGradient -->
 *   <exactLogCtmcRateGradient .../>
 * </benchmarkGradient>
 * }</pre>
 *
 * Plug this into a {@code <report>} block to trigger execution at the point
 * where BEAST normally prints model reports.
 *
 * @author Filippo Monti
 */
public class GradientBenchmarkParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "benchmarkGradient";
    public static final String WARMUP      = "warmup";
    public static final String ITERATIONS  = "iterations";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final int warmup     = xo.getAttribute(WARMUP,     50);
        final int iterations = xo.getAttribute(ITERATIONS, 100);

        final GradientWrtParameterProvider gradient =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        if (gradient == null) {
            throw new XMLParseException(
                    PARSER_NAME + ": must contain a GradientWrtParameterProvider child element");
        }

        return new GradientBenchmark(gradient, warmup, iterations);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(WARMUP,      true),
                AttributeRule.newIntegerRule(ITERATIONS,  true),
                new ElementRule(GradientWrtParameterProvider.class),
        };
    }

    @Override
    public String getParserDescription() {
        return "Benchmarks repeated gradient evaluations and reports mean/std/min/max per-call time.";
    }

    @Override
    public Class getReturnType() {
        return GradientBenchmark.class;
    }
}
