package dr.app.beauti.generator;

import dr.app.beauti.XMLWriter;
import dr.app.beauti.components.ComponentFactory;
import dr.app.beauti.options.BeautiOptions;
import dr.app.beauti.options.PartitionModel;
import dr.app.beauti.options.TreePrior;
import dr.evolution.util.Taxon;
import dr.evolution.util.TaxonList;
import dr.evomodel.speciation.SpeciationLikelihood;
import dr.evomodel.speciation.SpeciesBindings;
import dr.evomodel.speciation.SpeciesTreeBMPrior;
import dr.evomodel.speciation.SpeciesTreeModel;
import dr.evomodel.speciation.TreePartitionCoalescent;
import dr.evomodel.tree.TreeModel;
import dr.evomodelxml.BirthDeathModelParser;
import dr.evoxml.TaxonParser;
import dr.inference.distribution.ExponentialDistributionModel;
import dr.inference.model.ParameterParser;
import dr.util.Attribute;
import dr.xml.AttributeParser;
import dr.xml.XMLParser;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Alexei Drummond
 * @author Walter Xie
 */
public class MultiSpeciesCoalescentGenerator extends Generator {

    public MultiSpeciesCoalescentGenerator(BeautiOptions options, ComponentFactory[] components) {
        super(options, components);
    }
    
    /**
     * write tag <sp>
     * @param taxonList
     * @param writer
     */
    public void writeMultiSpecies(TaxonList taxonList, XMLWriter writer) {
    	List<String> species = options.getSpeciesList();
    	String sp;
    	
    	for (String eachSp : species) {
    		writer.writeOpenTag(SpeciesBindings.SP, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, eachSp)});

    		for (int i = 0; i < taxonList.getTaxonCount(); i++) {
    			Taxon taxon = taxonList.getTaxon(i);
    			sp = taxon.getAttribute(options.TRAIT_SPECIES).toString();

    			if (sp.equals(eachSp)) {
    				writer.writeIDref(TaxonParser.TAXON, taxon.getId());
    			}

    		}
    		writer.writeCloseTag(SpeciesBindings.SP);
    	}

    	writeGeneTrees (writer);
    }
    
    /**
     * write the species tree, species tree model, likelihood, etc.
     * @param writer
     */
    public void writeMultiSpeciesCoalescent(XMLWriter writer) {
    	writeSpeciesTree(writer);
    	writeSpeciesTreeModel(writer);
    	writeSpeciesTreeLikelihood(writer);
    	writeGeneUnderSpecies(writer);
    }
    
    
    private void writeGeneTrees(XMLWriter writer) {
    	writer.writeComment("Collection of Gene Trees");

        writer.writeOpenTag(SpeciesBindings.GENE_TREES, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, SpeciesBindings.GENE_TREES)});

    	// generate gene trees regarding each data partition
    	for (PartitionModel partitionModel : options.getActivePartitionModels()) {
    		writer.writeIDref(TreeModel.TREE_MODEL, partitionModel.getName() + "." + TreeModel.TREE_MODEL);
        }

        writer.writeCloseTag(SpeciesBindings.GENE_TREES);
    }


    private void writeSpeciesTree(XMLWriter writer) {
    	writer.writeComment("Species Tree: Provides Per branch demographic function");

        writer.writeOpenTag(SpeciesTreeModel.SPECIES_TREE, new Attribute[]{new Attribute.Default<String>(XMLParser.ID, SP_TREE),
				new Attribute.Default<String>(SpeciesTreeModel.BMPRIOR, "true")});        
        writer.writeIDref(options.TRAIT_SPECIES, options.TRAIT_SPECIES);
        
        //TODO: take sppSplitPopulations value from partionModel(?).constant.popSize
        // hard code get(0)
        double popSizeValue = options.getParameter("constant.popSize", options.getActivePartitionModels().get(0)).initial; // "initial" is "value"
        writer.writeOpenTag(SpeciesTreeModel.SPP_SPLIT_POPULATIONS, new Attribute[]{
        		new Attribute.Default<String>(AttributeParser.VALUE, Double.toString(popSizeValue))});
        
        writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
        		new Attribute.Default<String>(XMLParser.ID, SpeciesTreeModel.SPECIES_TREE + "." + SPLIT_POPS)}, true);
        
        writer.writeCloseTag(SpeciesTreeModel.SPP_SPLIT_POPULATIONS);
        
        writer.writeCloseTag(SpeciesTreeModel.SPECIES_TREE);

    }
    
    private void writeSpeciesTreeModel(XMLWriter writer) {    	    	
    	writer.writeComment("Species Tree: tree prior");
    	
    	if (options.nodeHeightPrior == TreePrior.SPECIES_BIRTH_DEATH) {
    		writer.writeComment("Birth Death Model");
    		
	    	writer.writeOpenTag(BirthDeathModelParser.BIRTH_DEATH_MODEL, new Attribute[]{
	    			new Attribute.Default<String>(XMLParser.ID, BirthDeathModelParser.BIRTH_DEATH),
					new Attribute.Default<String>(XMLParser.Utils.UNITS, XMLParser.Utils.SUBSTITUTIONS)});      
	    	
	    	writer.writeOpenTag(BirthDeathModelParser.BIRTHDIFF_RATE);
	    	
	    	writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
	        		new Attribute.Default<String>(XMLParser.ID, options.TRAIT_SPECIES + "." + BirthDeathModelParser.BIRTHDIFF_RATE_PARAM_NAME),
	        		new Attribute.Default<String>(AttributeParser.VALUE, "1"),
	        		new Attribute.Default<String>(ParameterParser.LOWER, "0"),
	        		new Attribute.Default<String>(ParameterParser.UPPER, "1000000")}, true);
	    	
	    	writer.writeCloseTag(BirthDeathModelParser.BIRTHDIFF_RATE);
	    	    	
	    	writer.writeOpenTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);
	    	
	    	writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
	        		new Attribute.Default<String>(XMLParser.ID, options.TRAIT_SPECIES + "." + BirthDeathModelParser.RELATIVE_DEATH_RATE_PARAM_NAME),
	        		new Attribute.Default<String>(AttributeParser.VALUE, "0.5"),
	        		new Attribute.Default<String>(ParameterParser.LOWER, "0"),
	        		new Attribute.Default<String>(ParameterParser.UPPER, "1")}, true);
	    	
	    	writer.writeCloseTag(BirthDeathModelParser.RELATIVE_DEATH_RATE);
	    	    	
	    	writer.writeCloseTag(BirthDeathModelParser.BIRTH_DEATH_MODEL); 
//    	} else if (options.nodeHeightPrior == TreePrior.SPECIES_YULE) {
    		//TODO: YULE model.
    	}
    	
    }
    
    
    private void writeSpeciesTreeLikelihood(XMLWriter writer) {
    	writer.writeComment("Species Tree: Likelihood of species tree");
	    	
    	if (options.nodeHeightPrior == TreePrior.SPECIES_BIRTH_DEATH) {
    		writer.writeComment("Birth Death Model");
    		
	    	writer.writeOpenTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD, new Attribute[]{
	    			new Attribute.Default<String>(XMLParser.ID, SPECIATION_LIKE)});      
	    	
	    	writer.writeOpenTag(SpeciationLikelihood.MODEL); 
	    	writer.writeIDref(BirthDeathModelParser.BIRTH_DEATH_MODEL, BirthDeathModelParser.BIRTH_DEATH);    	
	    	writer.writeCloseTag(SpeciationLikelihood.MODEL); 
	    	
	    	writer.writeOpenTag(SpeciesTreeModel.SPECIES_TREE); 
	    	writer.writeIDref(SpeciesTreeModel.SPECIES_TREE, SP_TREE);    	
	    	writer.writeCloseTag(SpeciesTreeModel.SPECIES_TREE); 
	    	
	    	writer.writeCloseTag(SpeciationLikelihood.SPECIATION_LIKELIHOOD); 
//    	} else if (options.nodeHeightPrior == TreePrior.SPECIES_YULE) {
    		//TODO: YULE model.
    	}
    }
    
    private void writeGeneUnderSpecies(XMLWriter writer) {
    	
    	writer.writeComment("Species Tree: Coalescent likelihood for gene trees under species tree");
    	
    	// speciesCoalescent id="coalescent"
    	writer.writeOpenTag(TreePartitionCoalescent.SPECIES_COALESCENT, new Attribute[]{
    			new Attribute.Default<String>(XMLParser.ID, COALESCENT)});   
    	
    	writer.writeIDref(options.TRAIT_SPECIES, options.TRAIT_SPECIES); 
    	writer.writeIDref(SpeciesTreeModel.SPECIES_TREE, SP_TREE); 
    	
    	writer.writeCloseTag(TreePartitionCoalescent.SPECIES_COALESCENT); 
    	
    	// exponentialDistributionModel id="pdist"
    	writer.writeOpenTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, new Attribute[]{
    			new Attribute.Default<String>(XMLParser.ID, PDIST)});  
    	
    	writer.writeOpenTag(ExponentialDistributionModel.MEAN); 
    	
    	writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
        		new Attribute.Default<String>(XMLParser.ID, options.TRAIT_SPECIES + "." + options.POP_MEAN),
        		new Attribute.Default<String>(AttributeParser.VALUE, "0.001")}, true);
    	
    	writer.writeCloseTag(ExponentialDistributionModel.MEAN); 
    	
    	writer.writeCloseTag(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL); 
    	
    	// STPopulationPrior id="stp" log_root="true"
    	writer.writeOpenTag(SpeciesTreeBMPrior.STPRIOR, new Attribute[]{
    			new Attribute.Default<String>(XMLParser.ID, STP),
    			new Attribute.Default<String>(SpeciesTreeBMPrior.LOG_ROOT, "true")});  
    	writer.writeIDref(SpeciesTreeModel.SPECIES_TREE, SP_TREE); 
    	
    	writer.writeOpenTag(SpeciesTreeBMPrior.TIPS);     	
    	
    	writer.writeIDref(ExponentialDistributionModel.EXPONENTIAL_DISTRIBUTION_MODEL, PDIST);
    	
    	writer.writeCloseTag(SpeciesTreeBMPrior.TIPS); 
    	
    	writer.writeOpenTag(SpeciesTreeBMPrior.STSIGMA);     	
    	
    	writer.writeTag(ParameterParser.PARAMETER, new Attribute[]{
    			// <parameter id="stsigma" value="1" /> 
        		new Attribute.Default<String>(XMLParser.ID, SpeciesTreeBMPrior.STSIGMA.toLowerCase()),
        		new Attribute.Default<String>(AttributeParser.VALUE, "1")}, true);
    	
    	writer.writeCloseTag(SpeciesTreeBMPrior.STSIGMA); 
    	
    	writer.writeCloseTag(SpeciesTreeBMPrior.STPRIOR); 
    }
}
