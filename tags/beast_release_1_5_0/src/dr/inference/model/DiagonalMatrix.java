package dr.inference.model;

import dr.xml.*;

/**
 * @author Marc Suchard
 */
public class DiagonalMatrix extends MatrixParameter {

	public final static String MATRIX_PARAMETER = "diagonalMatrix";

	private Parameter diagonalParameter;

	public DiagonalMatrix(Parameter param) {
		super(MATRIX_PARAMETER);
		addParameter(param);
		diagonalParameter = param;
	}

//	public DiagonalMatrix(String name, Parameter parameter) {
//		Parameter.Default(name, parameters);
//	}

	public double getParameterValue(int row, int col) {
		if (row != col)
			return 0.0;
		return diagonalParameter.getParameterValue(row);
	}

	public double[][] getParameterAsMatrix() {
		final int I = getDimension();
		double[][] parameterAsMatrix = new double[I][I];
		for (int i = 0; i < I; i++) {
			parameterAsMatrix[i][i] = diagonalParameter.getParameterValue(i);
		}
		return parameterAsMatrix;
	}

	public int getColumnDimension() {
		return diagonalParameter.getDimension();
	}

	public int getRowDimension() {
		return diagonalParameter.getDimension();
	}

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MATRIX_PARAMETER;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {


			Parameter diagonalParameter = (Parameter) xo.getChild(Parameter.class);

			return new DiagonalMatrix(diagonalParameter);
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A diagonal matrix parameter constructed from its diagonals.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(Parameter.class, 1, 1),
		};

		public Class getReturnType() {
			return MatrixParameter.class;
		}
	};


}
