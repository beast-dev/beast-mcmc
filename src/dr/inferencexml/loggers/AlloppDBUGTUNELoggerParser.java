package dr.inferencexml.loggers;

import java.io.PrintWriter;

import dr.evomodel.speciation.AlloppSpeciesBindings;
import dr.evomodel.speciation.AlloppSpeciesNetworkModel;
import dr.inference.loggers.AlloppDBUGTUNELogger;
import dr.util.FileHelpers;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.ElementRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;


/**
 * 
 * @author Graham Jones
 *         Date: 01/07/2011
 */


public class AlloppDBUGTUNELoggerParser extends AbstractXMLObjectParser {
	public static final String ALLOPPDBUGTUNE = "logAlloppDBUGTUNE";
	public static final String FILE_NAME = FileHelpers.FILE_NAME;
	
	@Override
	public String getParserName() {
		return ALLOPPDBUGTUNE;
	}

	@Override
	public Object parseXMLObject(XMLObject xo) throws XMLParseException {
		final PrintWriter pw = LoggerParser.getLogFile(xo, getParserName());
		final AlloppSpeciesNetworkModel apspnet = (AlloppSpeciesNetworkModel) xo.getChild(AlloppSpeciesNetworkModel.class);
		final AlloppSpeciesBindings apspb = (AlloppSpeciesBindings) xo.getChild(AlloppSpeciesBindings.class);
		return new AlloppDBUGTUNELogger(pw, apspnet, apspb);
	}

	@Override
	public XMLSyntaxRule[] getSyntaxRules() {
        return new XMLSyntaxRule[]{
        		new StringAttributeRule(FILE_NAME, "File path for debug info"),
                new ElementRule(AlloppSpeciesBindings.class),
                new ElementRule(AlloppSpeciesNetworkModel.class),
        };
	}

	@Override
	public String getParserDescription() {
		return "Logger for debugging and tuning information";
	}

	@Override
	public Class getReturnType() {
		return AlloppDBUGTUNELoggerParser.class;
	}

}
