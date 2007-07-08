package dr.inference.model;

import dr.xml.*;
//import org.w3c.dom.Document;
//import org.w3c.dom.Element;
//import dr.xml.*;

/**
 * Created by IntelliJ IDEA.
 * User: msuchard
 * Date: Dec 29, 2006
 * Time: 11:50:21 AM
 * To change this template use File | Settings | File Templates.
 */
public class MatrixParameter extends CompoundParameter {

	public final static String MATRIX_PARAMETER = "matrixParameter";

	public MatrixParameter(String name) {
		super(name);
	}

	public MatrixParameter(String name, Parameter[] parameters) {
		super(name, parameters);
	}

	public double getParameterValue(int row, int col) {
		return getParameter(col).getParameterValue(row);
	}

	public double[][] getParameterAsMatrix() {
		final int I = getRowDimension();
		final int J = getColumnDimension();
		double[][] parameterAsMatrix = new double[I][J];
		for (int i = 0; i < I; i++) {
			for (int j = 0; j < J; j++)
				parameterAsMatrix[i][j] = getParameterValue(i, j);
		}
		return parameterAsMatrix;
	}

	public int getColumnDimension() {
		return getNumberOfParameters();
	}

	public int getRowDimension() {
		return getParameter(0).getDimension();
	}

	// **************************************************************
	// XMLElement IMPLEMENTATION
	// **************************************************************

//    public Element createElement(Document d) {
//        throw new RuntimeException("Not implemented yet!");
//    }

	public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {

		public String getParserName() {
			return MATRIX_PARAMETER;
		}

		public Object parseXMLObject(XMLObject xo) throws XMLParseException {

			MatrixParameter matrixParameter = new MatrixParameter(MATRIX_PARAMETER);

			int dim = 0;

			for (int i = 0; i < xo.getChildCount(); i++) {
				Parameter parameter = (Parameter) xo.getChild(i);
				matrixParameter.addParameter(parameter);
				if (i == 0)
					dim = parameter.getDimension();
				else if (dim != parameter.getDimension())
					throw new XMLParseException("All parameters must have the same dimension to construct a rectangular matrix");
			}
			return matrixParameter;
		}

		//************************************************************************
		// AbstractXMLObjectParser implementation
		//************************************************************************

		public String getParserDescription() {
			return "A matrix parameter constructed from its component parameters.";
		}

		public XMLSyntaxRule[] getSyntaxRules() {
			return rules;
		}

		private XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
				new ElementRule(Parameter.class, 1, Integer.MAX_VALUE),
		};

		public Class getReturnType() {
			return MatrixParameter.class;
		}
	};


}
