package dr.inferencexml.distribution;

import dr.inference.distribution.IndexedLogitNormalWorkingPrior;
import dr.inference.model.IndexedParameter;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 * 
 * Example XML
 *
 * <indexedWorkingLogitNormalPriorLikelihood id="workingPrior">
 *   <ctsParameter>
 *     <parameter idref="ctsParameter"/>
 *   </ctsParameter>
 *   <atomsValues>
 *     <parameter idref="alphaAtoms"/>
 *   </atomsValues>
 *   <indicator>
 *     <parameter idref="indicator"/>
 *   </indicator>
 *   <atomsMapping>
 *     <indexedParameter idref="atomsMapping"/>
 *   </atomsMapping>
 * </indexedWorkingLogitNormalPriorLikelihood>
 */
public class IndexedLogitNormalWorkingPriorParser extends AbstractXMLObjectParser {

    public static final String NAME = IndexedLogitNormalWorkingPrior.ID;

    private static final String CTS_PARAMETER = "ctsParameter";
    private static final String ATOMS_VALUES = "atomsValues";
    private static final String INDICATOR  = "indicator";
    private static final String ATOM_INDEX = "atomsMapping";

    @Override
    public String getParserName() {
        return NAME;
    }

    @Override
    public Object parseXMLObject(final XMLObject xo) throws XMLParseException {

        final Parameter ctsParameter = (Parameter) xo.getChild(CTS_PARAMETER).getChild(Parameter.class);
        final Parameter atomsValues = (Parameter) xo.getChild(ATOMS_VALUES).getChild(Parameter.class);
        final Parameter indicator = (Parameter) xo.getChild(INDICATOR).getChild(Parameter.class);
        final IndexedParameter atomsMapping = (IndexedParameter) xo.getChild(ATOM_INDEX).getChild(IndexedParameter.class);

        if (indicator.getDimension() != ctsParameter.getDimension()) {
            throw new XMLParseException("<indicator> must have same dimension as <ctsParameter>.");
        }
        if (atomsMapping.getDimension() != ctsParameter.getDimension()) {
            throw new XMLParseException("<atomsMapping> must have same dimension as <ctsParameter>.");
        }
        
        return new IndexedLogitNormalWorkingPrior(ctsParameter, atomsValues, indicator, atomsMapping);
    }

    @Override
    public String getParserDescription() {
        return "Working prior likelihood for rho_b in (0,1) via logit-normal: "
                + "if z_b=0 then logit(rho_b) ~ Normal(0,1); "
                + "if z_b=1 then logit(rho_b) ~ Normal(logit(alpha[s_b]),1), with Jacobian included.";
    }

    @Override
    public Class getReturnType() {
        return IndexedLogitNormalWorkingPrior.class;
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
            new ElementRule(CTS_PARAMETER, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, false),
            new ElementRule(ATOMS_VALUES, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, false),
            new ElementRule(INDICATOR, new XMLSyntaxRule[] { new ElementRule(Parameter.class) }, false),
            new ElementRule(ATOM_INDEX, new XMLSyntaxRule[] { new ElementRule(IndexedParameter.class) }, false),
    };
}