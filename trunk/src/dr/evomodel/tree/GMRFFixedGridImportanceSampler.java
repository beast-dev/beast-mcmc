package dr.evomodel.tree;

import java.io.File;
import java.io.FileReader;
import java.io.Reader;

import dr.xml.AbstractXMLObjectParser;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class GMRFFixedGridImportanceSampler {

	public static final String FIXED_GRID_IMPORTANCE_SAMPLER = "gmrfFixedGridImportanceSampler";
	public static final String PROPOSAL_FILE_NAME = "proposals";
	public static final String LIKELIHOOD_FILE_NAME = "likelihood";
	public static final String BURN_IN = "burnIn";
	
	private Reader[] reader;
	private int burnIn;
	
	public GMRFFixedGridImportanceSampler(Reader[] r, int bI){
		reader = r;
		burnIn = bI;
	}
	
	public void report(){
		//Do something
	}
	
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "Generates importance samples from a fixed grid";
		}

		public Class getReturnType() {
			return GMRFFixedGridImportanceSampler.class;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
			
		};
		
		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			Reader[] reader;
			String proposalFileName = xo.getStringAttribute(PROPOSAL_FILE_NAME);
			String likelihoodFileName = xo.getStringAttribute(LIKELIHOOD_FILE_NAME);
			
			try{
				File proposalFile = new File(proposalFileName);
				
				String proposalName = proposalFile.getName();
				String proposalParent = proposalFile.getParent();

				if (!proposalFile.isAbsolute()) {
					proposalParent = System.getProperty("user.dir");
				}
				
				File likelihoodFile = new File(likelihoodFileName);
				
				String likelihoodName = likelihoodFile.getName();
				String likelihoodParent = likelihoodFile.getParent();

				if (!likelihoodFile.isAbsolute()) {
					likelihoodParent = System.getProperty("user.dir");
				}
				reader = new Reader[2];
				
				reader[0] = new FileReader(new File(proposalParent,proposalName));
				reader[1] = new FileReader(new File(likelihoodParent,likelihoodName));
								
			}catch(java.io.IOException e){
				//throw new XMLParseException("File '" + fileName + "' can not be opened for " + getParserName() + " element.");
				throw new XMLParseException("Something is wrong");
			}
			
			int burnin = -1;
			if (xo.hasAttribute(BURN_IN)) {
				burnin = xo.getIntegerAttribute(BURN_IN);
			}
			
			GMRFFixedGridImportanceSampler analysis = new GMRFFixedGridImportanceSampler(reader,burnin);
			
			analysis.report();
			
			System.out.println("");
			System.out.flush();
			
			return analysis;
		}

		public String getParserName() {
			// TODO Auto-generated method stub
			return null;
		}
		
	};
}
