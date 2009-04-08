package dr.evomodel.coalescent;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.util.Date;
import java.util.logging.Logger;

import dr.app.beast.BeastVersion;
import dr.evomodelxml.LoggerParser;
import dr.inference.loggers.LogColumn;
import dr.inference.loggers.LogFormatter;
import dr.inference.loggers.MCLogger;
import dr.inference.loggers.TabDelimitedFormatter;
import dr.math.MathUtils;
import dr.math.distributions.LogNormalDistribution;
import dr.xml.AbstractXMLObjectParser;
import dr.xml.AttributeRule;
import dr.xml.ElementRule;
import dr.xml.XMLObject;
import dr.xml.XMLObjectParser;
import dr.xml.XMLParseException;
import dr.xml.XMLSyntaxRule;

public class BayesianSkylineFixedGridLogger extends AbstractFixedGridLogger{

	public static final String BAYESIAN_SKYLINE_FIXED_GRID_LOGGER = "skylineFixedGridLogger";
	
	private BayesianSkylineLikelihood bsl;
	
	private int extraDraws = 10;
	private double extraDrawStdDev = 2.0;
	private double extraDensity = 0;
		
	public BayesianSkylineFixedGridLogger(LogFormatter formatter, int logEvery,
			BayesianSkylineLikelihood bsl, double gridHeight, int intervalNumber) {
		super(formatter, logEvery, gridHeight, intervalNumber);
				
		this.bsl = bsl;
	}
		
	public double[] getPopSizes(){
		
		extraDensity = 0;
		
		double[] bslPopSizes = bsl.getPopSizeParameter().getParameterValues();
		
		double[] popSizes = new double[bslPopSizes.length + extraDraws];
		
		for(int i = 0; i < bslPopSizes.length; i++){
			popSizes[i] = bslPopSizes[i];
		}
				
		for(int i = bslPopSizes.length; i < popSizes.length; i++){
			double extraDraw = MathUtils.nextGaussian()*extraDrawStdDev + Math.log(popSizes[i-1]);
						
			popSizes[i] = Math.exp(extraDraw);
			
			extraDensity += LogNormalDistribution.logPdf(popSizes[i], Math.log(popSizes[i-1]), extraDrawStdDev);
		}
				
		return popSizes;
	}
	
	public double[] getCoalescentHeights(){
		double[] bslGroupHeights = bsl.getGroupHeights();
		
		double[] groupHeights = new double[bslGroupHeights.length + extraDraws];
		
		for(int i = 0; i < bslGroupHeights.length; i++){
			groupHeights[i] = bslGroupHeights[i];
		}
		
		double length = (getGridStopTime() - bslGroupHeights[bslGroupHeights.length - 1])/extraDraws;
		
		for(int i = bslGroupHeights.length; i < groupHeights.length; i++){
			groupHeights[i] = groupHeights[i - 1] + length;
		}
		
		groupHeights[groupHeights.length - 1] = getGridStopTime();
		
		return  groupHeights;
	}
	
	public double getAdditionalDensity(){
		return extraDensity;
	}
		
	public static XMLObjectParser PARSER = new AbstractXMLObjectParser(){

		public String getParserDescription() {
			return "";
		}

		public Class getReturnType() {
			return BayesianSkylineFixedGridLogger.class;
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return new XMLSyntaxRule[]{
					new ElementRule(BayesianSkylineLikelihood.class),
					AttributeRule.newDoubleRule(GRID_STOP_TIME,false),
					AttributeRule.newIntegerRule(NUMBER_OF_INTERVALS,false),
					AttributeRule.newIntegerRule(LoggerParser.LOG_EVERY,true),
					AttributeRule.newStringRule(LoggerParser.FILE_NAME, true),
			};
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {
			
			int logEvery = 1;
			if(xo.hasAttribute(LoggerParser.LOG_EVERY)){
				logEvery = xo.getIntegerAttribute(LoggerParser.LOG_EVERY);
			}
			
			String fileName = null;
			if (xo.hasAttribute(LoggerParser.FILE_NAME)) {
				fileName = xo.getStringAttribute(LoggerParser.FILE_NAME);
			}
			
			double gridHeight = xo.getDoubleAttribute(GRID_STOP_TIME);
			int intervalNumber = xo.getIntegerAttribute(NUMBER_OF_INTERVALS);
						
			PrintWriter pw;
			
			Logger.getLogger("dr.evomodel").info("Creating " + BAYESIAN_SKYLINE_FIXED_GRID_LOGGER + " \n"
					+ "\t" + GRID_STOP_TIME + ": " + gridHeight  
					+ "\n\t" + NUMBER_OF_INTERVALS + ": " + intervalNumber);
			
			if (fileName != null) {

				try {
					File file = new File(fileName);
					String name = file.getName();
					String parent = file.getParent();

					if (!file.isAbsolute()) {
						parent = System.getProperty("user.dir");
					}
					pw = new PrintWriter(new FileOutputStream(new File(parent, name)));
				} catch (FileNotFoundException fnfe) {
					throw new XMLParseException("File '" + 
							fileName + "' can not be opened for " + 
							getParserName() + " element.");
				}
			} else {
				pw = new PrintWriter(System.out);
			}
			
			LogFormatter lf = new TabDelimitedFormatter(pw);
			
			MCLogger logger = new BayesianSkylineFixedGridLogger(lf,logEvery,
					(BayesianSkylineLikelihood)xo.getChild(BayesianSkylineLikelihood.class),
					gridHeight, intervalNumber);

			final BeastVersion version = new BeastVersion();
			
			logger.setTitle("BEAST " + version.getVersionString() + ", " + 
								version.getBuildString() + "\n" + 
								"Generated " + (new Date()).toString() + 
								"\nString = " + (MathUtils.getSeed()) +
								"\nFirst value corresponds to coalescent interval closet to sampling time\n" +
								"Last value corresponds to coalescent interval closet to the root\n" + 
								"Grid Height = " + gridHeight +  
								"\nNumber of intervals = " + intervalNumber);
								
			for(int i = 0 ; i < intervalNumber; i++){
				logger.addColumn(new LogColumn.Default("V" + (i+1),null));
			}
			logger.addColumn(new LogColumn.Default("AdditionalDensity",null));
			
			return (BayesianSkylineFixedGridLogger)logger;
			
		}

		public String getParserName() {
			return BAYESIAN_SKYLINE_FIXED_GRID_LOGGER;
		}
		
	};

}
