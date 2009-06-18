package dr.evomodel.coalescent;

import dr.inference.model.Statistic;
import dr.xml.*;

/**
 * @author Joseph Heled
 * @version $Id$
 *          Created by IntelliJ IDEA.
 *          User: joseph Date: 5/02/2007 Time: 16:17:11
 */
public class PopulationSizeGraph extends Statistic.Abstract {

    public static String POPGRAPH_STATISTIC = "popGraph";

    private double tm = 0;
    private VariableDemographicModel vdm = null;

    public PopulationSizeGraph(VariableDemographicModel vdm, double  tm) {
        super("popGraph");
        this.vdm = vdm;
        this.tm = tm;
    }

    public int getDimension() { return 1; }

    public double getStatisticValue(int dim) {
        return vdm.getDemographicFunction().getDemographic(tm);
    }

    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() { return POPGRAPH_STATISTIC; }

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

            Object child = xo.getChild(0);

            if (child instanceof VariableDemographicModel ) {
                final double dim = xo.getDoubleAttribute("time");
                return  new PopulationSizeGraph((VariableDemographicModel)child, dim);
            }

            throw new XMLParseException("Unknown element found in " + getParserName() + " element:" + child);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "This element returns a statistic that is the population size at evenly spaced intervals over tree.";
		}

		public Class getReturnType() { return PopulationSizeGraph.class; }

		public XMLSyntaxRule[] getSyntaxRules() { return rules; }

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
			new ElementRule(VariableDemographicModel.class, 1, 1 ),
                AttributeRule.newDoubleRule("time")
        };
	};
}
