package dr.evomodel.continuous;

import dr.inference.model.Parameter;
import dr.math.MathUtils;
import dr.xml.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;


public class MapDiffusionModel extends MultivariateDiffusionModel {
	public static final String MAP_DIFFUSION_MODEL = "mapDiffusionModel";
	public static final String GRASS_FILE_NAME = "grassMapFile";
	public static final String STARTING_VALUES = "randomStartingValues";

	public MapDiffusionModel(TopographicalMap map, Parameter graphRate) {
		super();
		this.map = map;
		this.graphRate = graphRate;
		addVariable(graphRate);
		initializationReport();
	}


	public TopographicalMap getMap() {
		return map;
	}

	private void initializationReport() {

		System.out.println("Constructing map diffusion model for");
		System.out.println("\tMap: " + map.getXDim() + "x" + map.getYDim());
		System.out.println("\tRandom-walk order: " + map.getOrder());
		System.out.println("\tRandom-walk non-zero size: " + map.getNonZeroSize());
//		System.out.println("\tRW1 Matrix:\n"+map.getMatrix().sparseRepresentation());
		System.out.println("\tRate parameter: " + graphRate.getStatisticName());
		//System.out.println("\tTest probability: "+map.getMatrix().getExponentialEntry(0,0,1.0));

	}

	/**
	 * @return the log likelihood of going from start to stop in the given time
	 */
	public double getLogLikelihood(double[] start, double[] stop, double time) {
//		System.err.println("Calc from "+start[0]+" "+start[1]);
//		System.err.println("Calc to   "+stop[0]+" "+stop[1]);

//		System.err.println("like");
		double rtn =

//		return
				Math.log(
						map.getCTMCProbability(start, stop, graphRate.getParameterValue(0) * time)
				);

//		System.err.println("logProb = "+rtn);
		return rtn;

	}


	public void handleParameterChangedEvent(Parameter parameter, int index) {
//		calculatePrecisionInfo();
	}

//	private void calculatePrecisionInfo() {
//	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MAP_DIFFUSION_MODEL;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {


			Parameter graphRate = (Parameter) xo.getChild(Parameter.class);

			String mapFileName = xo.getStringAttribute(GRASS_FILE_NAME);

			File mapFile;

			try {
				File file = new File(mapFileName);
				String name = file.getName();
				String parent = file.getParent();

				if (!file.isAbsolute()) {
					parent = System.getProperty("user.dir");
				}
				mapFile = new File(parent, name);
				new FileReader(mapFile);
			} catch (FileNotFoundException fnfe) {
				throw new XMLParseException("File '" + mapFileName + "' can not be opened for " + getParserName() + " element.");
			}

			double[][] mapValues;

			try {
				mapValues = TopographicalMap.readGRASSAscii(mapFileName);
			} catch (IOException e) {
				throw new XMLParseException("File '" + mapFileName + "' can not be read as GRASS file");
			}

			TopographicalMap map = new TopographicalMap(mapValues);

			XMLObject cxo = xo.getChild(STARTING_VALUES);
			if (cxo != null) {
				System.err.println("Init");
				Parameter startPosition = (Parameter) cxo.getChild(Parameter.class);
				int numberPositions = startPosition.getDimension() / 2;
				int dimX = map.getXDim();
				int dimY = map.getYDim();
				for (int i = 0; i < numberPositions; i++) {
					int x = -1;
					int y = -1;
					while (map.getIndex(x, y) == -1) {
						x = MathUtils.nextInt(dimX);
						y = MathUtils.nextInt(dimY);
						x = 0;
						y = 100;
					}
					startPosition.setParameterValue(i * 2, x);

					startPosition.setParameterValue(i * 2 + 1, y);
					System.err.println("set: " + x + "," + y);
				}
			}

			return new MapDiffusionModel(map, graphRate);

		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "Describes a multivariate discrete diffusion process.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private final XMLSyntaxRule[] rules = {
				AttributeRule.newStringRule(GRASS_FILE_NAME),
				new ElementRule(Parameter.class),
		};

		public Class getReturnType() {
			return MapDiffusionModel.class;
		}

	};

	private final TopographicalMap map;
	private final Parameter graphRate;

}
