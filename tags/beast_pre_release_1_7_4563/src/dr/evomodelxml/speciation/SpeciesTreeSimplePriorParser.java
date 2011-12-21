package dr.evomodelxml.speciation;

import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodel.speciation.SpeciesTreeSimplePrior;
import dr.inference.distribution.ParametricDistributionModel;
import dr.inference.model.Parameter;
import dr.xml.*;

/**
 */
public class SpeciesTreeSimplePriorParser extends AbstractXMLObjectParser {
    private static final String STPRIOR = "speciesTreePopulationPrior";
    public static final String TIPS = "tipsDistribution";

        public String getParserDescription() {
            return "";
        }

        public Class getReturnType() {
            return SpeciesTreeSimplePrior.class;
        }

        public String getParserName() {
            return STPRIOR;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
            SpeciesTreeModel st = (SpeciesTreeModel) xo.getChild(SpeciesTreeModel.class);

            //ParametricDistributionModel pr = (ParametricDistributionModel) xo.getChild(ParametricDistributionModel.class);
            Parameter pr = (Parameter)((XMLObject)xo.getChild("sigma")).getChild(Parameter.class);

            final XMLObject cxo = xo.getChild(TIPS);
            final ParametricDistributionModel tipsPrior = (ParametricDistributionModel) cxo.getChild(ParametricDistributionModel.class);

            return new SpeciesTreeSimplePrior(st, pr, tipsPrior);
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return new XMLSyntaxRule[]{
                    new ElementRule(SpeciesTreeModel.class),
                    new ElementRule(TIPS,
                            new XMLSyntaxRule[] { new ElementRule(ParametricDistributionModel.class) }),
                    //new ElementRule(ParametricDistributionModel.class),
                    new ElementRule("sigma", new XMLSyntaxRule[] { new ElementRule(Parameter.class) })
            };
        }
}
