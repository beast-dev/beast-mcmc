package dr.evomodel.bigfasttree.thorney;

import java.util.ArrayList;


public interface MutationList{
    public double getMutationCount();



    public class DetailedMutationList implements MutationList{
    
        private ArrayList<Mutation> mutations;

        public DetailedMutationList(){
            this.mutations = new ArrayList<>();
        }

        public double  getMutationCount(){
            return (double) mutations.size();
        };

        public int getSite(int mutationIndex){
            return mutations.get(mutationIndex).site;
        }

        public int getRef(int mutationIndex){
            return mutations.get(mutationIndex).ref;
        }

        public int getAlt(int mutationIndex){
            return mutations.get(mutationIndex).alt;
        }

        public void addMutation(Mutation mut){
            this.mutations.add(mut);
        }

        public Mutation removeMutation(int mutationIndex){
            return this.mutations.remove(mutationIndex);
        }


        protected class Mutation{

            final private int alt;
            final private int ref;
            final private int site;

            protected Mutation(int site, int ref, int alt){
                this.alt= alt;
                this.ref = ref;
                this.site = site;
            }
        }
    }

    public class SimpleMutationList implements MutationList{
        private double mutations;
        public SimpleMutationList(double muts){
            this.mutations = muts;
        }

        public double  getMutationCount(){
            return this.mutations;
        };
        public void setMutationCount(double muts){
            this.mutations=muts;
        }
    }
}
