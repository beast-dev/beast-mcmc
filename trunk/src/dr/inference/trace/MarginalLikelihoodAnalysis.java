package dr.inference.trace;

import dr.xml.*;
import dr.util.Attribute;
import dr.util.NumberFormatter;

import java.io.Reader;
import java.io.File;
import java.io.FileReader;
import java.io.FileNotFoundException;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Jan 17, 2007
 * Time: 6:45:12 PM
 * To change this template use File | Settings | File Templates.
 */
public class MarginalLikelihoodAnalysis {

	public static final String ML_ANALYSIS = "marginalLikelihoodAnalysis";
	public static final String FILE_NAME = "fileName";
	public static final String BURN_IN = "burnIn";
	public static final String COLUMN_NAME = "likelihoodColumn";

	private Trace trace;
	private int burnin;

	public MarginalLikelihoodAnalysis(Trace trace, int burnin) {
		this.trace = trace;
		this.burnin = burnin;
	}

	public double calculateLogMarginalLikelihood() {
		double sample[] = trace.getValues(burnin);
	}



	public String getParserName() { return ML_ANALYSIS; }

	public Object parseXMLObject(XMLObject xo) throws XMLParseException {

		try {
            Reader reader = null;

            String fileName = xo.getStringAttribute(FILE_NAME);
            try {
                File file = new File(fileName);
                String name = file.getName();
                String parent = file.getParent();

                if (!file.isAbsolute()) {
                    parent = System.getProperty("user.dir");
                }

                reader = new FileReader(new File(parent, name));
            } catch (FileNotFoundException fnfe) {
                throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
            }

			int burnin = -1;
			if (xo.hasAttribute(BURN_IN)) {
				// leaving the burnin attribute off will result in 10% being used
				burnin = xo.getIntegerAttribute(BURN_IN);
			}

//			TraceAnalysis[] analysis = TraceAnalysis.analyzeLogFile(reader,burnin);

//			Trace[] traces = Trace.Utils.loadTraces(reader);

			XMLObject child = (XMLObject)xo.getChild(COLUMN_NAME);
			String likelihoodName = child.getStringAttribute(Attribute.NAME);
			Trace trace = Trace.Utils.loadTrace(reader,likelihoodName);
			reader.close();

			int maxState = Trace.Utils.getMaximumState(trace);

			if (burnin == -1) {
				burnin = maxState / 10;
			} /*else if (percentage) {
				burnin = (int)Math.round((double)maxState * (double)burnin / 100.0);
			}*/

			if (burnin < 0 || burnin >= maxState) {
				burnin = maxState / 10;
				System.out.println("WARNING: Burn-in larger than total number of states - using to 10%");
			}

//					report(reader, burnin);
=
			MarginalLikelihoodAnalysis analysis = new  MarginalLikelihoodAnalysis(trace,burnin);



		} catch (java.io.IOException ioe) {
			throw new XMLParseException(ioe.getMessage());
		}
	}

	//************************************************************************
	// AbstractXMLObjectParser implementation
	//************************************************************************

	public String getParserDescription() {
		return "Performs a trace analysis. Estimates the mean of the various statistics in the given log file.";
	}

	public Class getReturnType() { return TraceAnalysis[].class; }

	public XMLSyntaxRule[] getSyntaxRules() { return rules; }

	private XMLSyntaxRule[] rules = new XMLSyntaxRule[] {
		new StringAttributeRule(FILE_NAME, "The name of a BEAST log file (can not include trees, which should be logged separately" ),
		AttributeRule.newIntegerRule("burnIn", true)
			//, "The number of states (not sampled states, but actual states) that are discarded from the beginning of the trace before doing the analysis" ),
	};
}
