package dr.evomodel.bigfasttree.thorney;


import dr.evolution.datatype.DataType;
import dr.evolution.tree.NodeRef;

public interface MutationBranchMap{
     public DataType getDataType();
     public double getMutations(final NodeRef node);

     
     public abstract class AbstractMutationBranchMap  implements MutationBranchMap{
          private DataType dataType;

     public  AbstractMutationBranchMap(DataType dataType){
          this.dataType = dataType;
     }

     public DataType getDataType() {
          return dataType;
     }
}
}
