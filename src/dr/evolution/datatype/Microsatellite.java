/*
 * Microsatellite.java
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

package dr.evolution.datatype;
import java.util.ArrayList;

/**
 * @author Chieh-Hsi Wu
 *
 * Microsatellite data type
 */
public class Microsatellite extends DataType {

    public static final String DESCRIPTION = "microsatellite";
    public static int UNKNOWN_STATE_LENGTH = -1;
    private int min;
    private int max;
    private int unitLength;
    private String name;
    public static final Microsatellite INSTANCE = new Microsatellite();

    public Microsatellite() {}

    /**
     * Constructor
     *
     * @param min   integer representing the minimum length of the microsatellite
     * @param max   integer representing the maximum length of the microsatellite
     */
    public Microsatellite(String name, int min, int max){
        this(name, min, max, 1);
    }

    public Microsatellite(int min, int max){
        this("microsat", min, max, 1);
    }

    /**
     * Constructor
     *
     * @param min           integer representing the minimum length of the microsatellite
     * @param max           integer representing the maximum length of the microsatellite
     * @param unitLength    the length in nucleotide units for one repeat unit.
     */
    public Microsatellite(String name, int min, int max, int unitLength){
        this.min = min;
        this.max = max;
        this.name = name;
        this.unitLength = unitLength;
        stateCount = (max - min)/unitLength + 1;
        ambiguousStateCount = stateCount + 1;
    }
     /**
      * This class constructs a Microsatellite Parameter.
      *
      * @param pattern      an array of ints representing the lengths (in nucleotide units) at a microsatellite locus
      * @param extRange     the length (in repeat units) to which to extend the current range of microsatellite in the pattern
      * @param unitLength   the length in nucleotide units for one repeat unit.
      */
    public Microsatellite(int[] pattern, int extRange, int unitLength){

         this.unitLength = unitLength;

         //find the current max and min in the pattern parameter
         min = pattern[0];
         max = pattern[0];
         for(int i = 1; i < pattern.length; i++){
             if(min > pattern[i]&&(pattern[i] > 0)){
                 min=pattern[i];
             }
             if(max < pattern[i]&&(pattern[i] > 0)){
                 max=pattern[i];
             }
         }

         //new max and min according to the extRange provided
         max = max + extRange*this.unitLength;
         min = min - extRange*this.unitLength;

         //set stateCount
         if((max-min)%this.unitLength == 0){
             stateCount = (max - min)/this.unitLength + 1;
         }else{
             throw new IllegalArgumentException("Incorrect microsatellite unit length.");
         }

         //set ambiguous stateCount
         ambiguousStateCount = stateCount + 1;

    }

    /**
     *
     * @param srtRawLength  a String representing the raw length of a microsatellite allele in nucleotide units
     * @return int          the state of microsatellite allele corresponding to the length
     */
    public int getState(String srtRawLength){
        char chRawLength =  srtRawLength.charAt(0);
        try{
            if(chRawLength == UNKNOWN_CHARACTER){
                return getState(UNKNOWN_STATE_LENGTH);
            }else{
                return getState(Integer.parseInt(srtRawLength));
            }
        } catch(java.lang.NumberFormatException exp) {
            throw new java.lang.NumberFormatException(srtRawLength+" can not be converted. State needs to be an integer or unknown (?).");
        }

    }

    /**
     *
     * @param rawLength an int representing the raw length of a microsatellite allele in nucleotide units
     * @return int      the state of microsatellite allele corresponding to the length
     */
    public int getState(int rawLength){
        if(rawLength != UNKNOWN_STATE_LENGTH){
            if (rawLength < min) {
                throw new java.lang.IllegalArgumentException("Microsatellite length value is less, (" + rawLength + ") than the specified minimum (" + min + ").");
            }
            return (int)Math.ceil(((double)rawLength - min)/unitLength);
        }else{
            return stateCount;
        }
    }


    /**
     *
     * @param strStates an ArrayList of Strings representing the raw microsatellite lengths in a pattern at a locus
     * @return int[]    an array of ints representing the raw microsatellite lengths in a pattern at a locus
     */
    public static int[] convertToLengths(ArrayList<String> strStates){
        return convertToLengths(strStates.toArray(new String[strStates.size()]));
    }

    /**
     *
     * @param strStates an array of Strings representing the raw microsatellite lengths in a pattern at a locus
     * @return int[]    an array of ints representing the raw microsatellite lengths in a pattern at a locus
     */
    public static int[] convertToLengths(String[] strStates){

        int[] lengths = new int[strStates.length];
        for(int i = 0; i < strStates.length; i++){
            char ch = strStates[i].charAt(0);
            try{
                if(ch == UNKNOWN_CHARACTER){
                    lengths[i] = UNKNOWN_STATE_LENGTH;
                }else{
                    lengths[i] = Integer.parseInt(strStates[i]);
                }
            }catch(java.lang.NumberFormatException exp){
                throw new java.lang.NumberFormatException(strStates[i]+" can not be converted. State needs to be an integer or unknown (?).");
            }
        }
        return lengths;
    }

    /**
     *
     * @param   stateCode   the code of the state of a microsatellite allele
     * @return  the raw length of the microsallite allele given it's state code
     *
     */
    public int getActualLength(int stateCode){
        if(stateCode <  stateCount){
            return (stateCode+min);
        }else if(stateCode == stateCount){
            return UNKNOWN_STATE_LENGTH;
        }else{
            throw new java.lang.RuntimeException("The given state must be an integer greater or equal to -1");
        }
    }

    /**
     *
     * @param  state    code representing a microsatellite allele
     * @return true if min <= state <= max, false otherwise.
     *
     */
    public boolean isWithinRange(int state){
        return (state >= min && state <= max);
    }

    /**
     * @return the upper bound of  allele length
     */
    public int getMax(){
        return max;
    }

    /**
     *
     * @return the lower bound of allele length
     *
     */
    public int getMin(){
        return min;
    }

    public void setMax(int max){
        this.max = max;
    }

    public void setMin(int min){
        this.min = min;
    }

    /**
     * @return true if the provided stateCode = unknownStateCode
     */
    public boolean isUnknownState(int stateCode){
        return (stateCode == stateCount);
    }

    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
     * @return number of unambiguous states for this data type.
     */
    public int getStateCount(){
        return stateCount;
    }

    /**
     * @return number of states (ambiguous states inclusive) for this data type.
     */
    public int getAmbiguousStateCount(){
        return ambiguousStateCount;
    }

    /**
     * @return the length (in bp) of one repeat unit
     */
    public int getUnitLength(){
        return unitLength;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getName() {
        return name;
    }    

    /**
     * @return the description of the data type
     */
    public String getDescription() {
		return DESCRIPTION;
	}

    public int getType(){
        return MICRO_SAT;
    }

}
