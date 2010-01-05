package dr.evomodel.arg;

import dr.inference.model.Statistic;
import dr.xml.*;


/**
 * @author Marc Suchard
 */
public class ARGReassortmentNodeCountStatistic extends Statistic.Abstract {

    public static final String REASSORTMENT_STATISTIC = "argReassortmentNodeCount";


    public ARGReassortmentNodeCountStatistic(String name, ARGModel arg) {
        super(name);
        this.arg = arg;      
    }

    public int getDimension() {
        return 1;
    }

    public double getStatisticValue(int dim) {
        return arg.getReassortmentNodeCount();
    }

     public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

        public String getParserName() {
            return REASSORTMENT_STATISTIC;
        }

        public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            String name = xo.getAttribute(NAME, xo.getId());
            ARGModel arg = (ARGModel) xo.getChild(ARGModel.class);

            return new ARGReassortmentNodeCountStatistic(name,arg);
        }

        //************************************************************************
        // AbstractXMLObjectParser implementation
        //************************************************************************

        public String getParserDescription() {
            return "A statistic that returns the number of reassortment nodes in an ARG";
        }

        public Class getReturnType() {
            return Statistic.class;
        }

        public XMLSyntaxRule[] getSyntaxRules() {
            return rules;
        }

        private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
                new StringAttributeRule(NAME, "A name for this statistic for the purpose of logging", true),
                new ElementRule(ARGModel.class),
        };

    };


    private ARGModel arg;

}
