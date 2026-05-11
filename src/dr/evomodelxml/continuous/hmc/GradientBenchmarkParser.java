package dr.evomodelxml.continuous.hmc;

import dr.evomodel.treedatalikelihood.discrete.GradientBenchmark;
import dr.inference.hmc.GradientWrtParameterProvider;
import dr.inference.model.Likelihood;
import dr.xml.*;

/**
 * XML parser for {@link GradientBenchmark}.
 *
 * Syntax:
 * <pre>{@code
 * <benchmarkGradient warmup="500" iterations="200" dirty="true">
 *   <exactLogCtmcRateGradient idref="spectralExactGradient"/>
 *   <treeDataLikelihood idref="treeLikelihood"/>
 * </benchmarkGradient>
 * }</pre>
 *
 * When {@code dirty="true"}, {@code makeDirty()} is called on the provided
 * {@code Likelihood} before each iteration, simulating the full HMC per-step
 * cost (pre-order + post-order traversal + Fréchet kernel).
 *
 * @author Filippo Monti
 */
public class GradientBenchmarkParser extends AbstractXMLObjectParser {

    public static final String PARSER_NAME = "benchmarkGradient";
    public static final String WARMUP      = "warmup";
    public static final String ITERATIONS  = "iterations";
    public static final String DIRTY       = "dirty";

    @Override
    public String getParserName() { return PARSER_NAME; }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        final int warmup     = xo.getAttribute(WARMUP,     50);
        final int iterations = xo.getAttribute(ITERATIONS, 100);
        final boolean dirty  = xo.getAttribute(DIRTY, false);

        final GradientWrtParameterProvider gradient =
                (GradientWrtParameterProvider) xo.getChild(GradientWrtParameterProvider.class);
        if (gradient == null) {
            throw new XMLParseException(
                    PARSER_NAME + ": must contain a GradientWrtParameterProvider child element");
        }

        final Likelihood dirtyTarget = (Likelihood) xo.getChild(Likelihood.class);
        if (dirty && dirtyTarget == null) {
            throw new XMLParseException(
                    PARSER_NAME + ": dirty=\"true\" requires a Likelihood child element (e.g. treeDataLikelihood)");
        }

        return new GradientBenchmark(gradient, warmup, iterations, dirty, dirtyTarget);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
                AttributeRule.newIntegerRule(WARMUP,      true),
                AttributeRule.newIntegerRule(ITERATIONS,  true),
                AttributeRule.newBooleanRule(DIRTY,       true),
                new ElementRule(GradientWrtParameterProvider.class),
                new ElementRule(Likelihood.class, true),
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
