/*
 * CollectionHash.java
 *
 * Copyright (C) 2002-2006 Alexei Drummond and Andrew Rambaut
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

package dr.util;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Enumeration;
import java.util.Hashtable;


// todo not used?
//
//public class CollectionHash {
//
//    private final Hashtable table;
//
//    public CollectionHash() { table = new Hashtable(); }
//
//    public void put(Object key, Object o) {
//		Collection c = (Collection)table.get(key);
//
//		if (c != null) {
//		    c.add(o);
//		} else {
//		    Collection newc = new ArrayList();
//		    newc.add(o);
//		    table.put(key, newc);
//		}
//    }
//
//    public Object get(Object key) {
//		Collection c = (Collection)table.get(key);
//
//		if (c == null) return null;
//
//		return c.iterator().next();
//    }
//
//    public Enumeration keys() {
//		return table.keys();
//    }
//
//    public Object[] getAll(Object key) {
//		return getCollection(key).toArray();
//    }
//
//    public Collection getCollection(Object key) {
//		return (Collection)table.get(key);
//    }
//
//    public int getSize(Object key) {
//        Collection collection = (Collection)table.get(key);
//        if (collection == null) return 0;
//        return collection.size();
//    }
//
//    public void remove(Object key) {
//    	table.remove(key);
//    }
//}
