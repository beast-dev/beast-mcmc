package dr.evolution;

import dr.util.Identifiable;

import java.util.ArrayList;

/**
 * @author Aaron Darling (koadman)
 */
public class LinkageConstraints implements Identifiable {

	protected String id = null;
	ArrayList<LinkedGroup> linkedGroups;

	public LinkageConstraints(ArrayList<LinkedGroup> linkedGroups){
		this.linkedGroups = linkedGroups;
	}

	public ArrayList<LinkedGroup> getLinkedGroups(){
		return linkedGroups;
	}
	

	public String getId() {
		return id;
	}

	
	public void setId(String id) {
		this.id = id;
	}

}
