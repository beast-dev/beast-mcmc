package dr.evolution;

import java.util.ArrayList;

import dr.util.Identifiable;

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
	
	@Override
	public String getId() {
		return id;
	}

	@Override
	public void setId(String id) {
		this.id = id;
	}

}
