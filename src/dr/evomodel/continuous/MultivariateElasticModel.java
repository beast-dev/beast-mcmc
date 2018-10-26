/*
 * MultivariateElasticModel.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package dr.evomodel.continuous;

import dr.evolution.tree.Tree;
import dr.evolution.tree.TreeAttributeProvider;
import dr.inference.model.*;
import dr.math.matrixAlgebra.missingData.MissingOps;
import org.ejml.data.Complex64F;
import org.ejml.data.DenseMatrix64F;
import org.ejml.factory.DecompositionFactory;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.ops.EigenOps;

/**
 * @author Marc Suchard
 * @author Paul Bastide
 */


public class MultivariateElasticModel extends AbstractModel implements TreeAttributeProvider {

    private static final String ELASTIC_PROCESS = "multivariateElasticModel";
    //    public static final String ELASTIC_CONSTANT = "strengthOfSelectionMatrix";
    private static final String ELASTIC_TREE_ATTRIBUTE = "strengthOfSelection";

    public static final double LOG2PI = Math.log(2 * Math.PI);

    /**
     * Construct a diffusion model.
     */

    public MultivariateElasticModel(MatrixParameterInterface strengthOfSelectionMatrixParameter) {

        super(ELASTIC_PROCESS);

        this.strengthOfSelectionMatrixParameter = strengthOfSelectionMatrixParameter;
        dim = strengthOfSelectionMatrixParameter.getRowDimension();
        assert dim == strengthOfSelectionMatrixParameter.getColumnDimension() : "Strength of Selection matrix should be square.";

        if (strengthOfSelectionMatrixParameter instanceof DiagonalMatrix) {
            this.parametrization = Parametrization.AS_DIAGONAL;
        } else if (strengthOfSelectionMatrixParameter instanceof CompoundEigenMatrix) {
            this.parametrization = Parametrization.AS_DECOMPOSED;
        } else {
            this.parametrization = Parametrization.GENERAL;
        }

        isSymmetric = strengthOfSelectionMatrixParameter.isConstrainedSymmetric();

        calculateSelectionInfo();
        addVariable(strengthOfSelectionMatrixParameter);

    }

    public MultivariateElasticModel() {
        super(ELASTIC_PROCESS);
    }

    private void calculateSelectionInfo() {
//        strengthOfSelectionMatrix = strengthOfSelectionMatrixParameter.getParameterAsMatrix();
        this.eigenDecompositionStrengthOfSelection = parametrization.decomposeStrenghtOfSelection(strengthOfSelectionMatrixParameter, dim, isSymmetric);
    }

    public MatrixParameterInterface getStrengthOfSelectionMatrixParameter() {
        checkVariableChanged();
        return strengthOfSelectionMatrixParameter;
    }

    public double[][] getStrengthOfSelectionMatrix() {
        if (strengthOfSelectionMatrixParameter != null) {
            checkVariableChanged();
            return strengthOfSelectionMatrixParameter.getParameterAsMatrix();
        }
        return null;
    }

    public double[] getStrengthOfSelectionMatrixAsVector() {
        double[][] strengthOfSelectionMatrix = getStrengthOfSelectionMatrix();
        int dim = strengthOfSelectionMatrix.length;
        double[] strengthOfSelectionVector = new double[dim * dim];
        for (int i = 0; i < dim; ++i) {
            System.arraycopy(strengthOfSelectionMatrix[i], 0, strengthOfSelectionVector, i * dim, dim);
        }
        return strengthOfSelectionVector;
    }

    public double[] getEigenValuesStrengthOfSelection() {
        if (strengthOfSelectionMatrixParameter != null) {
            checkVariableChanged();
            return parametrization.eigenValuesMatrix(strengthOfSelectionMatrixParameter, eigenDecompositionStrengthOfSelection, dim);
        }
        return null;
    }

    public double[] getEigenVectorsStrengthOfSelection() {
        if (strengthOfSelectionMatrixParameter != null) {
            checkVariableChanged();
            return parametrization.eigenVectorsMatrix(strengthOfSelectionMatrixParameter, eigenDecompositionStrengthOfSelection);
        }
        return null;
    }

    private void checkVariableChanged() {
        if (variableChanged) {
            calculateSelectionInfo();
            variableChanged = false;
        }
    }

    public boolean isDiagonal() {
        return parametrization == Parametrization.AS_DIAGONAL;
    }

    public boolean isSymmetric() {
        return isSymmetric;
    }

    // *****************************************************************
    // Parametrization
    // *****************************************************************
    enum Parametrization {
        AS_DIAGONAL {
            @Override
            public EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface AParam, int dim, boolean isSymmetric) {
                return null;
            }
            @Override
            public double[] eigenValuesMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA, int dim) {
                return ((DiagonalMatrix) AParam).getDiagonalParameter().getParameterValues();
            }
            @Override
            public double[] eigenVectorsMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA) {
                return null;
            }
        },
        AS_DECOMPOSED {
            @Override
            public EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface AParam, int dim, boolean isSymmetric) {
                return null;
            }
            @Override
            public double[] eigenValuesMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA, int dim) {
                return ((CompoundEigenMatrix) AParam).getEigenValues();
            }
            @Override
            public double[] eigenVectorsMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA) {
                return ((CompoundEigenMatrix) AParam).getEigenVectors();
            }
        },
        GENERAL {
            @Override
            public EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface AParam, int dim, boolean isSymmetric) {
                DenseMatrix64F A = MissingOps.wrap(AParam);
                EigenDecomposition eigA = DecompositionFactory.eig(dim, true, isSymmetric);
                if (!eigA.decompose(A)) throw new RuntimeException("Eigen decomposition failed.");
                return eigA;
            }
            @Override
            public double[] eigenValuesMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA, int dim) {
                assert eigDecompA != null : "The eigen decomposition should already be computed at this point.";
                double[] eigA = new double[dim];
                for (int p = 0; p < dim; ++p) {
                    Complex64F ev = eigDecompA.getEigenvalue(p);
                    assert ev.isReal() : "Selection strength A should only have real eigenvalues.";
                    assert ev.real > 0 : "Selection strength A should only have positive real eigenvalues.";
                    eigA[p] = ev.real;
                }
                return eigA;
            }
            @Override
            public double[] eigenVectorsMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA) {
                assert eigDecompA != null : "The eigen decomposition should already be computed at this point.";
                return EigenOps.createMatrixV(eigDecompA).getData();
            }
        };

        abstract EigenDecomposition decomposeStrenghtOfSelection(MatrixParameterInterface AParam, int dim, boolean isSymmetric);
        abstract double[] eigenValuesMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA, int dim);
        abstract double[] eigenVectorsMatrix(MatrixParameterInterface AParam, EigenDecomposition eigDecompA);
    }

    // *****************************************************************
    // Interface Model
    // *****************************************************************
    public void handleModelChangedEvent(Model model, Object object, int index) {
        // no intermediates need to be recalculated...
    }

    protected void handleVariableChangedEvent(Variable variable, int index, Parameter.ChangeType type) {
        variableChanged = true;
    }

    protected void storeState() {
//        savedStrengthOfSelectionMatrix = strengthOfSelectionMatrix;
        savedEigenDecompositionStrengthOfSelection = eigenDecompositionStrengthOfSelection;
        storedVariableChanged = variableChanged;
    }

    protected void restoreState() {
//        strengthOfSelectionMatrix = savedStrengthOfSelectionMatrix;
        eigenDecompositionStrengthOfSelection = savedEigenDecompositionStrengthOfSelection;
        variableChanged = storedVariableChanged;
    }

    protected void acceptState() {
    } // no additional state needs accepting

    public String[] getTreeAttributeLabel() {
        return new String[]{ELASTIC_TREE_ATTRIBUTE};
    }

    public String[] getAttributeForTree(Tree tree) {
        if (strengthOfSelectionMatrixParameter != null) {
            return new String[]{strengthOfSelectionMatrixParameter.toString()};
        }

        strengthOfSelectionMatrixParameter.toString();
        return new String[]{"null"};
    }

    // **************************************************************
    // XMLElement IMPLEMENTATION
    // **************************************************************

//    public Element createElement(Document document) {
//        throw new RuntimeException("Not implemented!");
//    }

//    // **************************************************************
//    // XMLObjectParser
//    // **************************************************************
//
//    public static XMLObjectParser PARSER = new AbstractXMLObjectParser() {
//
//        public String getParserName() {
//            return ELASTIC_PROCESS;
//        }
//
//        public Object parseXMLObject(XMLObject xo) throws XMLParseException {
//
//            XMLObject cxo = xo.getChild(ELASTIC_CONSTANT);
//            MatrixParameterInterface elasticParam = (MatrixParameterInterface)
//                    cxo.getChild(MatrixParameterInterface.class);
//
//            return new MultivariateElasticModel(elasticParam);
//        }
//
//        //************************************************************************
//        // AbstractXMLObjectParser implementation
//        //************************************************************************
//
//        public String getParserDescription() {
//            return "Describes a multivariate elastic process.";
//        }
//
//        public XMLSyntaxRule[] getSyntaxRules() {
//            return rules;
//        }
//
//        private final XMLSyntaxRule[] rules = {
//                new ElementRule(ELASTIC_CONSTANT,
//                        new XMLSyntaxRule[]{new ElementRule(MatrixParameterInterface.class)}),
//        };
//
//        public Class getReturnType() {
//            return dr.evomodel.continuous.MultivariateDiffusionModel.class;
//        }
//    };

    // **************************************************************
    // Private instance variables
    // **************************************************************

    private MatrixParameterInterface strengthOfSelectionMatrixParameter;
//    private double[][] strengthOfSelectionMatrix;
//    private double[][] savedStrengthOfSelectionMatrix;
    private EigenDecomposition eigenDecompositionStrengthOfSelection = null;
    private EigenDecomposition savedEigenDecompositionStrengthOfSelection = null;

    private Parametrization parametrization;
    private boolean isSymmetric;

    private int dim;

    private boolean variableChanged = true;
    private boolean storedVariableChanged;

}