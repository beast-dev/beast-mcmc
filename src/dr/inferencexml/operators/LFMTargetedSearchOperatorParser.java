package dr.inferencexml.operators;

import dr.inference.model.MatrixParameter;
import dr.inference.model.MatrixParameterInterface;
import dr.inference.model.Parameter;
import dr.inference.operators.LFMTargetedSearchOperator;
import dr.util.Transform;
import dr.xml.*;

import java.util.*;

public class LFMTargetedSearchOperatorParser extends AbstractXMLObjectParser {
    public static final String LFM_TARGETED_SEARCH_OPERATOR = "LFMTargetedSearchOperator";
    public static final String SPARSE_MATRIX = "sparseMatrix";
    public static final String LOADINGS_MATRIX = "loadingsMatrix";
    public static final String FACTORS_MATRIX = "factorsMatrix";
    public static final String SPARSE_TARGET_MATRICES = "sparseTargetMatrices";
    public static final String LOADINGS_TARGET_MATRICES = "loadingsTargetMatrices";
    public static final String FACTORS_TARGET_MATRICES = "factorsTargetMatrices";
    public static final String CUTOFFS = "cutoffs";
    public static final String WEIGHT = "weight";


    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        MatrixParameterInterface sparseMatrix = (MatrixParameterInterface) xo.getChild(SPARSE_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface loadingsMatrix = (MatrixParameterInterface) xo.getChild(LOADINGS_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface factorsMatrix = (MatrixParameterInterface) xo.getChild(FACTORS_MATRIX).getChild(MatrixParameterInterface.class);
        MatrixParameterInterface cutoffs = (MatrixParameterInterface) xo.getChild(CUTOFFS).getChild(MatrixParameterInterface.class);

        XMLObject cxo = xo.getChild(SPARSE_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> sparseTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < cxo.getChildCount(); i++) {
            sparseTargetList.add((MatrixParameterInterface) cxo.getChild(i));
        }
        XMLObject dxo = xo.getChild(LOADINGS_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> loadingsTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < dxo.getChildCount(); i++) {
            loadingsTargetList.add((MatrixParameterInterface) dxo.getChild(i));
        }
        XMLObject exo = xo.getChild(FACTORS_TARGET_MATRICES);
        ArrayList<MatrixParameterInterface> factorsTargetList = new ArrayList<MatrixParameterInterface>();
        for (int i = 0; i < exo.getChildCount(); i++) {
            factorsTargetList.add((MatrixParameterInterface) exo.getChild(i));
        }
        double weight = xo.getDoubleAttribute(WEIGHT);

        return new LFMTargetedSearchOperator(weight, sparseMatrix, sparseTargetList, factorsMatrix, factorsTargetList, loadingsMatrix, loadingsTargetList, cutoffs);
    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }
    XMLSyntaxRule[] rules = new XMLSyntaxRule[]{
            new ElementRule(SPARSE_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(SPARSE_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(FACTORS_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(FACTORS_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(LOADINGS_MATRIX, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
            new ElementRule(LOADINGS_TARGET_MATRICES, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class, 1, Integer.MAX_VALUE)}),
            new ElementRule(CUTOFFS, new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),

    };

    @Override
    public String getParserDescription() {
        return "Targeted search for Sparse Latent Factor Model";
    }

    @Override
    public Class getReturnType() {
        return LFMTargetedSearchOperator.class;
    }

    @Override
    public String getParserName() {
        return LFM_TARGETED_SEARCH_OPERATOR;
    }
}
