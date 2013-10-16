package dr.evomodelxml.speciation;

import dr.evolution.tree.Tree;
import dr.evomodel.speciation.MaskableSpeciationModel;
import dr.evomodel.speciation.ModelAveragingSpeciationLikelihood;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;
import dr.xml.*;

import java.util.ArrayList;
import java.util.List;
import java.util.logging.Logger;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class ModelAveragingSpeciationLikelihoodParser extends AbstractXMLObjectParser {
    public static final String MODEL_AVE_SPECIATION_LIKELIHOOD = "modelAveragingSpeciationLikelihood";
    public static final String MODEL = "sM";
    public static final String TREE = "sT";
    public static final String INDEX = "modelIndex";
    public static final String MAX_INDEX = "maxIndex";

    public String getParserName() {
        return MODEL_AVE_SPECIATION_LIKELIHOOD;
    }

    public Object parseXMLObject(XMLObject xo) throws XMLParseException {
        List<Tree> trees = new ArrayList<Tree>();
        List<MaskableSpeciationModel> models = new ArrayList<MaskableSpeciationModel>();
        Variable<Integer> index;

        System.out.println("id = " + xo.getId());

        XMLObject cxo = xo.getChild(MODEL);
        for (int m = 0; m < cxo.getChildCount(); m++) {
            final MaskableSpeciationModel specModel = (MaskableSpeciationModel) cxo.getChild(m);
            models.add(specModel);
        }

        cxo = xo.getChild(TREE);
        for (int t = 0; t < cxo.getChildCount(); t++) {
            final Tree tree = (Tree) cxo.getChild(t);
            trees.add(tree);
        }

//        cxo = xo.getChild(INDEX);
        index = (Variable<Integer>) xo.getElementFirstChild(INDEX); // integer index parameter size = real size - 1
        Parameter maxIndex = (Parameter) xo.getElementFirstChild(MAX_INDEX);
//        System.out.println(index.getClass());
//        for (int i=0; i<index.getSize(); i++) {
//            System.out.println(index.getValue(i).getClass());
//            index.setValue(i, 0);
//        }

        int indexLength = models.size();
        if (indexLength < 1 || trees.size() < 1) {
            throw new XMLParseException("It requires at least one tree or one speciation model.");
        } else if (indexLength != trees.size()) {
            throw new XMLParseException("The number of trees and the number of speciation models should be equal.");
        } else if (indexLength != index.getSize() + 1) { // integer index parameter size = real size - 1
            throw new XMLParseException("Index parameter must be same size as the number of trees.");
        }


        Logger.getLogger("dr.evomodel").info("Speciation model excluding " + " taxa remaining.");

        return new ModelAveragingSpeciationLikelihood(trees, models, index, maxIndex, xo.getId());
    }

    //************************************************************************
    // AbstractXMLObjectParser implementation
    //************************************************************************

    public String getParserDescription() {
        return "Model Averaging Speciation Likelihood.";
    }

    public Class getReturnType() {
        return ModelAveragingSpeciationLikelihood.class;
    }

    public XMLSyntaxRule[] getSyntaxRules() {
        return rules;
    }

    private final XMLSyntaxRule[] rules = {
            new ElementRule(MODEL, new XMLSyntaxRule[]{
                    new ElementRule(MaskableSpeciationModel.class, 1, Integer.MAX_VALUE)
            }),
            new ElementRule(TREE, new XMLSyntaxRule[]{
                    new ElementRule(Tree.class, 1, Integer.MAX_VALUE)
            }),
            new ElementRule(INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Variable.class)
            }),
            new ElementRule(MAX_INDEX, new XMLSyntaxRule[]{
                    new ElementRule(Parameter.class)
            }, true),
    };

}
