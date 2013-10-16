package dr.app.beagle.tools.parsers;

import dr.app.beagle.evomodel.treelikelihood.BeagleTreeLikelihood;
import dr.app.beagle.tools.SiteLogLikelihoodLogger;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class SiteLogLikelihoodLoggerParser extends AbstractXMLObjectParser {

	public static final String SITE_LOGLIKELIHOOD_LOGGER = "siteLogLikelihood";

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		SiteLogLikelihoodLogger siteLogLikelihoodLogger;
		BeagleTreeLikelihood beagleTreeLikelihood = null;

		for (int i = 0; i < xo.getChildCount(); i++) {
			beagleTreeLikelihood = (BeagleTreeLikelihood) xo.getChild(i);
		}

		siteLogLikelihoodLogger = new SiteLogLikelihoodLogger(
				beagleTreeLikelihood);

		return siteLogLikelihoodLogger;
	}// END: parseXMLObject

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
		return new XMLSyntaxRule[] { new ElementRule(BeagleTreeLikelihood.class) };
	}// END: getSyntaxRules

	@Override
	public String getParserName() {
		return SITE_LOGLIKELIHOOD_LOGGER;
	}// END: getParserName

	@Override
	public String getParserDescription() {
		return "Beagle site logLikelihood";
	}// END: getParserDescription

	@Override
	public Class<SiteLogLikelihoodLogger> getReturnType() {
		return SiteLogLikelihoodLogger.class;
	}// getReturnType

}// END: class
