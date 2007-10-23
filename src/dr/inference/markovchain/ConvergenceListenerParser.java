/**
 * 
 */
package dr.inference.markovchain;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.StringAttributeRule;
import dr.xml.XMLObject;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

/**
 * @author Sebastian Hoehna
 * 
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

      double threshold = 5.0;
      if (xo.hasAttribute(THRESHOLD)) {
         threshold = xo.getDoubleAttribute(THRESHOLD);
      }

      int checkEvery = 1;
      if (xo.hasAttribute(CHECK_EVERY)) {
         checkEvery = xo.getIntegerAttribute(CHECK_EVERY);
      }

      ConvergenceListener convergenceListener = new ConvergenceListener(null,
            checkEvery, threshold, treeFileName, referenceFileName);

      return convergenceListener;
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
