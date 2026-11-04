/*
 * MutationList.java
 *
 * Copyright © 2002-2024 the BEAST Development Team
 * http://beast.community/about
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
 *
 */

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
