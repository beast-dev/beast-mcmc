package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.AgeDependentBirthDeathIEModel;
import dr.inference.model.Parameter;
import dr.xml.*;

public class AgeDependentBirthDeathIEModelParser extends AbstractXMLObjectParser {

    private static final String PARSER_NAME = "ageDependentBirthDeathIEModel";
    private static final String ORIGIN_TIME = "originTime";
    private static final String EPOCH_TIMES = "epochTimes";
    private static final String BIRTH_SCALE = "birthScale";
    private static final String BIRTH_SHAPE = "birthShape";
    private static final String DEATH_SCALE = "deathScale";
    private static final String DEATH_SHAPE = "deathShape";

    private static final String STEPS = "timeSteps";
    private static final String EPS_PICARD = "epsPicard";
    private static final String MAX_ITER_PICARD = "maxIterPicard";
    private static final String SOLVER = "solver";
    private static final String THREADS = "threads";
    private static final String EXCLUDE_ROOT_BRANCH = "excludeRootBranch";

    public String getParserName() {
        return PARSER_NAME;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        Tree tree = (Tree) xo.getChild(Tree.class);

        double originTime = xo.getDoubleAttribute(ORIGIN_TIME);

        Parameter epochTimes = null;
        if (xo.hasChildNamed(EPOCH_TIMES)) {
            epochTimes = (Parameter) xo.getElementFirstChild(EPOCH_TIMES);
        }

        Parameter birthScale = (Parameter) xo.getElementFirstChild(BIRTH_SCALE);
        Parameter birthShape = (Parameter) xo.getElementFirstChild(BIRTH_SHAPE);
        Parameter deathScale = (Parameter) xo.getElementFirstChild(DEATH_SCALE);
        Parameter deathShape = (Parameter) xo.getElementFirstChild(DEATH_SHAPE);

        int numSteps = xo.getIntegerAttribute(STEPS);
        double epsPicard = xo.getDoubleAttribute(EPS_PICARD);
        int maxIterPicard = xo.getIntegerAttribute(MAX_ITER_PICARD);

        boolean useDirectQuadrature = false;
        if (xo.hasAttribute(SOLVER)) {
            String solver = xo.getStringAttribute(SOLVER);
            if ("quadrature".equalsIgnoreCase(solver)) {
                useDirectQuadrature = true;
            } else if (!"fft".equalsIgnoreCase(solver)) {
                throw new XMLParseException("solver must be 'fft' or 'quadrature', got: " + solver);
            }
        }

        int numThreads = xo.getAttribute(THREADS, 1);
        boolean excludeRootBranch = xo.getAttribute(EXCLUDE_ROOT_BRANCH, false);

        int numEpochs = (epochTimes != null) ? epochTimes.getDimension() + 1 : 1;

        if (birthScale.getDimension() != 1 && birthScale.getDimension() != numEpochs) {
            throw new XMLParseException("birthScale dimension (" + birthScale.getDimension() +
                    ") must be 1 or equal to epochTimes dimension (" + numEpochs + ")");
        }
        if (deathScale.getDimension() != 1 && deathScale.getDimension() != numEpochs) {
            throw new XMLParseException("deathScale dimension (" + deathScale.getDimension() +
                    ") must be 1 or equal to epochTimes dimension (" + numEpochs + ")");
        }
        if (birthShape.getDimension() != 2) {
            throw new XMLParseException("birthShape must have dimension 2 [b, gamma], got " + birthShape.getDimension());
        }
        if (deathShape.getDimension() != 2) {
            throw new XMLParseException("deathShape must have dimension 2 [b, gamma], got " + deathShape.getDimension());
        }

        return new AgeDependentBirthDeathIEModel(
                xo.getId(),
                tree,
                epochTimes,
                originTime,
                birthScale,
                birthShape,
                deathScale,
                deathShape,
                numSteps,
                epsPicard,
                maxIterPicard,
                useDirectQuadrature,
                numThreads,
                excludeRootBranch
        );
    }

    public String getParserDescription() {
        return "Optimized variant of ageDependentSkylineBirthDeathModel: same algorithm, " +
               "preallocated scratch pools, kernel caching, partial store/restore.";
    }

    public Class getReturnType() {
        return AgeDependentBirthDeathIEModel.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Tree.class),
            AttributeRule.newDoubleRule(ORIGIN_TIME),
            new ElementRule(EPOCH_TIMES, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }, true),
            new ElementRule(BIRTH_SCALE, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            new ElementRule(BIRTH_SHAPE, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            new ElementRule(DEATH_SCALE, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            new ElementRule(DEATH_SHAPE, new XMLSyntaxRule[]{ new ElementRule(Parameter.class) }),
            AttributeRule.newIntegerRule(STEPS),
            AttributeRule.newDoubleRule(EPS_PICARD),
            AttributeRule.newIntegerRule(MAX_ITER_PICARD),
            AttributeRule.newStringRule(SOLVER, true),
            AttributeRule.newIntegerRule(THREADS, true),
            AttributeRule.newBooleanRule(EXCLUDE_ROOT_BRANCH, true),
    };
}
