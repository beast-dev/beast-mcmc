/**
 *
 */
package dr.inference.markovchain;

import dr.xml.*;

/**
 * @author Sebastian Hoehna
 */
public class ConvergenceListenerParser extends AbstractXMLObjectParser {

    public final static String CONVERGENCE_LISTENER = "convergenceListener";

    public final static String TREE_FILE_NAME = "treeFilename";

    public final static String REFERENCE_FILE_NAME = "referenceFilename";

    public final static String THRESHOLD = "threshold";

    public final static String CHECK_EVERY = "checkEvery";

    /*
    * (non-Javadoc)
    * 
    * @see dr.xml.AbstractXMLObjectParser#getParserDescription()
    */
    @Override
    public String getParserDescription() {
        return "Checks the convergence of a run till it is under the threshold.";
    }

    /*
    * (non-Javadoc)
    * 
    * @see dr.xml.AbstractXMLObjectParser#getReturnType()
    */
    @Override
    public Class getReturnType() {
        return ConvergenceListener.class;
    }

    /*
    * (non-Javadoc)
    * 
    * @see dr.xml.AbstractXMLObjectParser#getSyntaxRules()
    */
    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new StringAttributeRule(TREE_FILE_NAME, "name of a tree log file", "trees.log"),
            new StringAttributeRule(REFERENCE_FILE_NAME, "name of a reference tree file", "trees.log"),
            AttributeRule.newIntegerRule(CHECK_EVERY, true),
            AttributeRule.newDoubleRule(THRESHOLD, true)
    };

    /*
    * (non-Javadoc)
    * 
    * @see dr.xml.AbstractXMLObjectParser#parseXMLObject(dr.xml.XMLObject)
    */
    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String treeFileName = xo.getStringAttribute(TREE_FILE_NAME);

        String referenceFileName = xo.getStringAttribute(REFERENCE_FILE_NAME);

        double threshold = xo.getAttribute(THRESHOLD, 5.0);

        int checkEvery = xo.getAttribute(CHECK_EVERY, 1);

        return new ConvergenceListener(null, checkEvery, threshold, treeFileName, referenceFileName);
    }

    /*
    * (non-Javadoc)
    * 
    * @see dr.xml.XMLObjectParser#getParserName()
    */
    public String getParserName() {
        return CONVERGENCE_LISTENER;
    }

}
