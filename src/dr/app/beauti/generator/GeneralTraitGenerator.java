package dr.app.beauti.generator;

import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.enumTypes.ClockType;
import dr.app.beauti.enumTypes.LocationSubstModelType;
import dr.app.beauti.options.*;
import dr.app.beauti.util.XMLWriter;
import dr.app.util.Arguments;
import dr.evolution.datatype.Nucleotides;
import dr.evolution.util.Taxon;
import dr.evomodel.branchratemodel.BranchRateModel;
import dr.evomodel.sitemodel.SiteModel;
import dr.evomodel.substmodel.AbstractSubstitutionModel;
import dr.evomodel.substmodel.SVSGeneralSubstitutionModel;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.branchratemodel.DiscretizedBranchRatesParser;
import dr.evomodelxml.branchratemodel.StrictClockBranchRatesParser;
import dr.evomodelxml.clock.ACLikelihoodParser;
import dr.evomodelxml.substmodel.ComplexSubstitutionModelParser;
import dr.evomodelxml.treelikelihood.AncestralStateTreeLikelihoodParser;
import dr.evomodelxml.treelikelihood.TreeLikelihoodParser;
import dr.evoxml.AttributePatternsParser;
import dr.evoxml.GeneralDataTypeParser;
import dr.evoxml.TaxaParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.util.List;

/**
 * @author Walter Xie
 */
public class GeneralTraitGenerator extends Generator {

    public static final String DATA = "Data";


    public GeneralTraitGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }

    public static String getLocationSubstModelTag(PartitionSubstitutionModel substModel) {
        if (substModel.getLocationSubstType() == LocationSubstModelType.SYM_SUBST) {
            return SVSGeneralSubstitutionModel.SVS_GENERAL_SUBSTITUTION_MODEL;
        } else if (substModel.getLocationSubstType() == LocationSubstModelType.ASYM_SUBST) {
            return ComplexSubstitutionModelParser.COMPLEX_SUBSTITUTION_MODEL;
        } else {
            return null;
        }
    }

    /**
     * write attr in taxon.
     *
     * @param taxon  Taxon
     * @param writer XMLWriter
     * @throws dr.app.util.Arguments.ArgumentException
     *          ArgumentException
     */
    public void writeAtrrTrait(Taxon taxon, XMLWriter writer) throws Arguments.ArgumentException {
        for (TraitData trait : BeautiOptions.getDiscreteIntegerTraits()) {
            if (!trait.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString())) {

                if (!taxon.containsAttribute(trait.getName())) {
                    throw new Arguments.ArgumentException("Cannot find trait " + trait.getName()
                            + "\nin taxon " + taxon.getId());
                }

                writer.writeOpenTag(AttributeParser.ATTRIBUTE, new Attribute[]{
                        new Attribute.Default<String>(Attribute.NAME, trait.getName())});

                writer.writeText(taxon.getAttribute(trait.getName()).toString());
                writer.writeCloseTag(AttributeParser.ATTRIBUTE);
            }
        }
    }

    /**
     * write <generalDataType> and <attributePatterns>
     *
     * @param traitData TraitData
     * @param writer    XMLWriter
     */
    public void writeGeneralDataType(TraitData traitData, XMLWriter writer) {
        writer.writeComment("trait = " + traitData.getName() + " trait_type = " + traitData.getTraitType());

        List<String> generalData = TraitData.getStatesListOfTrait(BeautiOptions.taxonList, traitData.getName());

        // <generalDataType>
        writer.writeOpenTag(GeneralDataTypeParser.GENERAL_DATA_TYPE, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitData.getPrefix() + DATA)});

        int numOfSates = generalData.size();
        writer.writeComment("Number Of Sates = " + numOfSates);

        for (String eachGD : generalData) {
            writer.writeTag(GeneralDataTypeParser.STATE, new Attribute[]{
                    new Attribute.Default<String>(GeneralDataTypeParser.CODE, eachGD)}, true);
        }

        writer.writeCloseTag(GeneralDataTypeParser.GENERAL_DATA_TYPE);

        // <attributePatterns>
        writer.writeOpenTag(AttributePatternsParser.ATTRIBUTE_PATTERNS, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitData.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS),
                new Attribute.Default<String>(AttributePatternsParser.ATTRIBUTE, traitData.getName())});
        writer.writeIDref(TaxaParser.TAXA, TaxaParser.TAXA);
        writer.writeIDref(GeneralDataTypeParser.GENERAL_DATA_TYPE, traitData.getPrefix() + DATA);
        writer.writeCloseTag(AttributePatternsParser.ATTRIBUTE_PATTERNS);
    }


    /**
     * Ancestral Tree Likelihood
     *
     * @param traitData TraitData
     * @param writer    XMLWriter
     */
    public void writeAncestralTreeLikelihood(TraitData traitData, XMLWriter writer) {
        PartitionSubstitutionModel substModel = traitData.getPartitionSubstitutionModel();
        PartitionTreeModel treeModel = traitData.getPartitionTreeModel();
        PartitionClockModel clockModel = traitData.getPartitionClockModel();

        writer.writeOpenTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD, new Attribute[]{
                new Attribute.Default<String>(XMLParser.ID, traitData.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD)});

        writer.writeIDref(AttributePatternsParser.ATTRIBUTE_PATTERNS, traitData.getPrefix() + AttributePatternsParser.ATTRIBUTE_PATTERNS);
        writer.writeIDref(TreeModel.TREE_MODEL, treeModel.getPrefix() + TreeModel.TREE_MODEL);
        writer.writeIDref(SiteModel.SITE_MODEL, substModel.getPrefix() + SiteModel.SITE_MODEL);
        writer.writeIDref(getLocationSubstModelTag(substModel), substModel.getPrefix() + AbstractSubstitutionModel.MODEL);

        switch (clockModel.getClockType()) {
            case STRICT_CLOCK:
                writer.writeIDref(StrictClockBranchRatesParser.STRICT_CLOCK_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;
            case UNCORRELATED_EXPONENTIAL:
            case UNCORRELATED_LOGNORMAL:
            case RANDOM_LOCAL_CLOCK:
                writer.writeIDref(DiscretizedBranchRatesParser.DISCRETIZED_BRANCH_RATES,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;

            case AUTOCORRELATED_LOGNORMAL:
                writer.writeIDref(ACLikelihoodParser.AC_LIKELIHOOD,
                        clockModel.getPrefix() + BranchRateModel.BRANCH_RATES);
                break;

            default:
                throw new IllegalArgumentException("Unknown clock model");
        }

        writer.writeCloseTag(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD);
    }

    public void writeAncestralTreeLikelihoodReferences(XMLWriter writer) {
        for (TraitData traitData : BeautiOptions.getDiscreteIntegerTraits()) { // Each TD except Species has one AncestralTreeLikelihood
            if (!traitData.getName().equalsIgnoreCase(TraitData.Traits.TRAIT_SPECIES.toString()))
                writer.writeIDref(AncestralStateTreeLikelihoodParser.RECONSTRUCTING_TREE_LIKELIHOOD,
                        traitData.getPrefix() + TreeLikelihoodParser.TREE_LIKELIHOOD);
        }
    }

    public void writeLogs(XMLWriter writer) {

    }

}

