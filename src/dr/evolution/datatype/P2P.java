/*
 * P2P.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
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
 */

package dr.evolution.datatype;

/**
 * @version 08/01/2010
 *
 * @author Marc A. Suchard
 * @author Yu-Nong Gong
 */
public class P2P extends DataType {

    /**
	 * Name of data type. For XML and human reading of data type.
	 */
	public static final String DESCRIPTION = "P2P";
	public static final int TYPE = 42;
	public static final P2P INSTANCE = new P2P();

//  public static final int ZERO_STATE = 0;
//	public static final int ONE_STATE = 1;

//	public static final int UNKNOWN_STATE = 2;
//	public static final int GAP_STATE = 3;

	/**
	 * A table to translate state numbers (0-3) into character codes
	 */
//	public static final char[] P2P_CHARS =
//		{ 'A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R',
//		'S','T','V','W','Y','B','Z','X'};

	/**
	 * A table to map state numbers (0-3) to their ambiguities
	 */
	public static final String[] P2P_AMBIGUITIES = {
"AA","AR","AN","AD","AC","AQ","AE","AG","AH","AI","AL","AK","AM","AF","AP","AS","AT","AW","AY","AV",
"RA","RR","RN","RD","RC","RQ","RE","RG","RH","RI","RL","RK","RM","RF","RP","RS","RT","RW","RY","RV",
"NA","NR","NN","ND","NC","NQ","NE","NG","NH","NI","NL","NK","NM","NF","NP","NS","NT","NW","NY","NV",
"DA","DR","DN","DD","DC","DQ","DE","DG","DH","DI","DL","DK","DM","DF","DP","DS","DT","DW","DY","DV",
"CA","CR","CN","CD","CC","CQ","CE","CG","CH","CI","CL","CK","CM","CF","CP","CS","CT","CW","CY","CV",
"QA","QR","QN","QD","QC","QQ","QE","QG","QH","QI","QL","QK","QM","QF","QP","QS","QT","QW","QY","QV",
"EA","ER","EN","ED","EC","EQ","EE","EG","EH","EI","EL","EK","EM","EF","EP","ES","ET","EW","EY","EV",
"GA","GR","GN","GD","GC","GQ","GE","GG","GH","GI","GL","GK","GM","GF","GP","GS","GT","GW","GY","GV",
"HA","HR","HN","HD","HC","HQ","HE","HG","HH","HI","HL","HK","HM","HF","HP","HS","HT","HW","HY","HV",
"IA","IR","IN","ID","IC","IQ","IE","IG","IH","II","IL","IK","IM","IF","IP","IS","IT","IW","IY","IV",
"LA","LR","LN","LD","LC","LQ","LE","LG","LH","LI","LL","LK","LM","LF","LP","LS","LT","LW","LY","LV",
"KA","KR","KN","KD","KC","KQ","KE","KG","KH","KI","KL","KK","KM","KF","KP","KS","KT","KW","KY","KV",
"MA","MR","MN","MD","MC","MQ","ME","MG","MH","MI","ML","MK","MM","MF","MP","MS","MT","MW","MY","MV",
"FA","FR","FN","FD","FC","FQ","FE","FG","FH","FI","FL","FK","FM","FF","FP","FS","FT","FW","FY","FV",
"PA","PR","PN","PD","PC","PQ","PE","PG","PH","PI","PL","PK","PM","PF","PP","PS","PT","PW","PY","PV",
"SA","SR","SN","SD","SC","SQ","SE","SG","SH","SI","SL","SK","SM","SF","SP","SS","ST","SW","SY","SV",
"TA","TR","TN","TD","TC","TQ","TE","TG","TH","TI","TL","TK","TM","TF","TP","TS","TT","TW","TY","TV",
"WA","WR","WN","WD","WC","WQ","WE","WG","WH","WI","WL","WK","WM","WF","WP","WS","WT","WW","WY","WV",
"YA","YR","YN","YD","YC","YQ","YE","YG","YH","YI","YL","YK","YM","YF","YP","YS","YT","YW","YY","YV",
"VA","VR","VN","VD","VC","VQ","VE","VG","VH","VI","VL","VK","VM","VF","VP","VS","VT","VW","VY","VV"
	};

	/**
	 * Private constructor - DEFAULT_INSTANCE provides the only instance
	 */
	private P2P() {
		stateCount = 400;
//		ambiguousStateCount = 4;
	}


    @Override
    public char[] getValidChars() {
        return null;
    }

    /**
	 * Get character corresponding to a given state
	 *
	 * @param state state
	 *
	 * return corresponding character
	 */
    public char getChar(int state) {
        throw new IllegalArgumentException("P2P datatype cannot be expressed as char");
    }
    
	public String getCode(int state) {
      return P2P_AMBIGUITIES[state];
	}

	/**
	 * returns an array containing the non-ambiguous states
	 * that this state represents.
	 */
//	public int[] getStates(int state) {
//		String stateString = P2P_AMBIGUITIES[state];
//        int[] states = new int[stateString.length()];
//		for (int i = 0; i < stateString.length(); i++) {
//			states[i] = getState(stateString.charAt(i));
//        }
//        return states;
//	}

	/**
	 * returns an array containing the non-ambiguous states that this state represents.
	 */
	public boolean[] getStateSet(int state) {

		boolean[] stateSet = new boolean[stateCount];
		for(int i=0;i<stateCount;i++){
			stateSet[i] = true;
		}

		return stateSet;
	}

	/**
	 * description of data type
	 *
	 * @return string describing the data type
	 */
	public String getDescription() {
		return DESCRIPTION;
	}

	/**
	 * type of data type
	 *
	 * @return integer code for the data type
	 */
	public int getType() {
		return TYPE;
	}

}
