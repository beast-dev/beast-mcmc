package dr.evolution;

import dr.evolution.util.TaxonList;
import dr.util.Identifiable;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkedGroup implements Identifiable {

	protected String id = null;
	TaxonList taxa;
	double probability;	// probability that reads in this group are linked to each other

	public LinkedGroup(TaxonList taxa, double probability){
		this.taxa = taxa;
		this.probability = probability;		
	}

	public LinkedGroup(TaxonList taxa){
		// default constructor assumes probability of linkage is 1
		this(taxa, 1.0);
	}
	
	public double getLinkageProbability(){
		return probability;
	}
	public TaxonList getLinkedReads(){
		return taxa;
	}


	public String getId() {
		return id;
	}

	
	public void setId(String id) {
		this.id = id;
	}

}
