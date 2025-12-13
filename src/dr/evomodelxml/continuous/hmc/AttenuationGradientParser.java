package dr.evomodelxml.continuous.hmc;

import dr.evolution.tree.Tree;
import dr.evomodel.treedatalikelihood.DataLikelihoodDelegate;
import dr.evomodel.treedatalikelihood.TreeDataLikelihood;
import dr.evomodel.treedatalikelihood.continuous.*;
import dr.evomodel.treedatalikelihood.continuous.backprop.*;
import dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient;
import dr.evomodelxml.treelikelihood.TreeTraitParserUtilities;
import dr.inference.model.*;
import dr.xml.*;

import java.util.ArrayList;
import java.util.Arrays;

import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDecomposedAttenuationGradient;
import static dr.evomodel.treedatalikelihood.hmc.AbstractDiffusionGradient.ParameterDiffusionGradient.createDiagonalAttenuationGradient;
import static dr.evomodelxml.treelikelihood.TreeTraitParserUtilities.DEFAULT_TRAIT_NAME;

/**
 * @author Paul Bastide
 * @author Marc A. Suchard
 */


public class AttenuationGradientParser extends AbstractXMLObjectParser {
    private final static String PRECISION_GRADIENT = "attenuationGradient";
    private final static String PARAMETER = "parameter";
    private final static String ATTENUATION_CORRELATION = "correlation";
    private final static String ATTENUATION_DIAGONAL = "diagonal";
    private final static String ATTENUATION_BOTH = "both";
    private static final String TRAIT_NAME = TreeTraitParserUtilities.TRAIT_NAME;
    private static final String APPROXIMATE = "approximate";


    @Override
    public String getParserName() {
        return PRECISION_GRADIENT;
    }

    private ParameterMode parseParameterMode(XMLObject xo) throws XMLParseException {
        // Choose which parameter(s) to update.
        ParameterMode mode = ParameterMode.WRT_BOTH;
        String parameterString = xo.getAttribute(PARAMETER, ATTENUATION_BOTH).toLowerCase();
        if (parameterString.compareTo(ATTENUATION_CORRELATION) == 0) {
            mode = ParameterMode.WRT_CORRELATION;
        } else if (parameterString.compareTo(ATTENUATION_DIAGONAL) == 0) {
            mode = ParameterMode.WRT_DIAGONAL;
        }
        return mode;
    }

    enum ParameterMode {
        WRT_BOTH {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter) {
                // Check if using BlockDiagonalCosSinMatrixParameter (backprop)
                if (parameter instanceof BlockDiagonalCosSinMatrixParameter) {
                    throw new RuntimeException("Gradient wrt full attenuation not yet implemented outside backpropagation.");
//                    return BackPropParameterDiffusionGradient.createBackpropAttenuationGradient(
//                            branchSpecificGradient,
//                            treeDataLikelihood,
//                            (BlockDiagonalCosSinMatrixParameter) parameter
//                    );
                }
                return createDecomposedAttenuationGradient(branchSpecificGradient, treeDataLikelihood, parameter);
//                throw new RuntimeException("Gradient wrt full attenuation not yet implemented.");
            }
        },
        WRT_CORRELATION {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter) {
                throw new RuntimeException("Gradient wrt correlation of attenuation not yet implemented.");
            }
        },
        WRT_DIAGONAL {
            @Override
            public Object factory(BranchSpecificGradient branchSpecificGradient,
                                  TreeDataLikelihood treeDataLikelihood,
                                  MatrixParameterInterface parameter) {
                return createDiagonalAttenuationGradient(branchSpecificGradient, treeDataLikelihood, parameter);
            }
        };

        abstract Object factory(BranchSpecificGradient branchSpecificGradient,
                                TreeDataLikelihood treeDataLikelihood,
                                MatrixParameterInterface parameter);
    }

    @Override
    public Object parseXMLObject(XMLObject xo) throws XMLParseException {

        String traitName = xo.getAttribute(TRAIT_NAME, DEFAULT_TRAIT_NAME);
        TreeDataLikelihood treeDataLikelihood = (TreeDataLikelihood) xo.getChild(TreeDataLikelihood.class);

        DataLikelihoodDelegate delegate = treeDataLikelihood.getDataLikelihoodDelegate();
        int dim = treeDataLikelihood.getDataLikelihoodDelegate().getTraitDim();
        Tree tree = treeDataLikelihood.getTree();

        ContinuousDataLikelihoodDelegate continuousData = (ContinuousDataLikelihoodDelegate) delegate;
        ContinuousTraitGradientForBranch traitGradient;

        CompoundEigenMatrix compoundEigenMatrix = (CompoundEigenMatrix) xo.getChild(CompoundEigenMatrix.class);
        BlockDiagonalCosSinMatrixParameter blockDiagonalCosSinMatrixParameter = (BlockDiagonalCosSinMatrixParameter) xo.getChild(BlockDiagonalCosSinMatrixParameter.class);
        if (compoundEigenMatrix != null) {
//            Parameter parameter = compoundEigenMatrix.getDiagonalParameter();
            traitGradient =
                    new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                            dim, tree, continuousData,
                            new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>(
                                    Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_GENERAL_SELECTION_STRENGTH)
                            ));
            BranchSpecificGradient branchSpecificGradient =
                    new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, compoundEigenMatrix);

            ParameterMode parameterMode = parseParameterMode(xo);
            return parameterMode.factory(branchSpecificGradient, treeDataLikelihood, compoundEigenMatrix);
        } else if (blockDiagonalCosSinMatrixParameter != null) {
            if (false) {
                traitGradient =
                        new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                                dim, tree, continuousData,
                                new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>(
                                        Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_GENERAL_SELECTION_STRENGTH)
                                ));
            } else {
                BlockDiagCosSinPrimitiveGradientMapper primitiveMapper =
                        new BlockDiagCosSinPrimitiveGradientMapper(blockDiagonalCosSinMatrixParameter);
                Parameter parameter = blockDiagonalCosSinMatrixParameter.getParameter();
                BackPropParameterProvider backPropParameterProvider = new BackPropParameterProvider.Default(treeDataLikelihood, continuousData, parameter, primitiveMapper);

                return backPropParameterProvider;
            }

            BranchSpecificGradient branchSpecificGradient =
                    new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, blockDiagonalCosSinMatrixParameter);

            ParameterMode parameterMode = parseParameterMode(xo);
            return parameterMode.factory(branchSpecificGradient, treeDataLikelihood, blockDiagonalCosSinMatrixParameter);
//            return new AttenuationBackpropGradient(branchSpecificGradient, treeDataLikelihood, blockDiagonalCosSinMatrixParameter, 10.0, 0.0);
        } else {

            MatrixParameterInterface parameter = (MatrixParameterInterface) xo.getChild(MatrixParameterInterface.class);
            if (parameter instanceof DiagonalMatrix) {
                traitGradient =
                        new ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient(
                                dim, tree, continuousData,
                                new ArrayList<ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter>(
                                        Arrays.asList(ContinuousTraitGradientForBranch.ContinuousProcessParameterGradient.DerivationParameter.WRT_DIAGONAL_SELECTION_STRENGTH)
                                ));
                BranchSpecificGradient branchSpecificGradient =
                        new BranchSpecificGradient(traitName, treeDataLikelihood, continuousData, traitGradient, parameter);

                ParameterMode parameterMode = parseParameterMode(xo);
                return parameterMode.factory(branchSpecificGradient, treeDataLikelihood, parameter);
            } else {
                throw new XMLParseException("Unsupported matrix parameter type for attenuation gradient: " + parameter.getClass().getName());
            }
        }

    }

    @Override
    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(Likelihood.class),
            new ElementRule(MatrixParameterInterface.class),
    };

    @Override
    public String getParserDescription() {
        return null;
    }

    @Override
    public Class getReturnType() {
        return AbstractDiffusionGradient.ParameterDiffusionGradient.class;
    }
}