/*
 * Vector3D.java
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

package dr.geo.math;

public final class Vector3D {

	/** Component X */
	protected double x;
	/** Component Y */
	protected double y;
	/** Component Z */
	protected double z;

	/** Origin Vector3D */
	public static final Vector3D ORIGIN = new Vector3D(0.0, 0.0, 0.0);

	/**
	* Create a Vector3D by components.
	* @param x component X
	* @param y component Y
	* @param z component Z
	*/
	public Vector3D(double x, double y, double z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}

	/**
	* Create a Vector3D from another Vector3D.
	*/
	public Vector3D(Vector3D a) {
		this.x = a.x;
		this.y = a.y;
		this.z = a.z;
	}
	    
	/**
	* add two Vector3Ds together.
	* @param a Vector3D to add to this one
	* @return a new <code>Vector3D</code>
	*/
	public Vector3D add(Vector3D a) {
		return (new Vector3D(x + a.x, y + a.y, z + a.z));
	}

	/**
	* Add two Vector3Ds (modifies this Vector3D).
	* @param a Vector3D to add
	* @return this <code>Vector3D</code>
	*/
	public Vector3D addU(Vector3D a) {
		x += a.x;
		y += a.y;
		z += a.z;
		return (this);
	}

	/**
	* Cross product of two Vector3Ds
	*/
	public Vector3D cross(Vector3D a) {
		return (new Vector3D(y * a.z - z * a.y, z * a.x - x * a.z, x * a.y - y * a.x));
	}

	/**
	* Dot product of two Vector3Ds
	*/
	public double dot(Vector3D a) {
		return (x * a.x + y * a.y + z * a.z);
	}
	
	public void negate() {
		x = -x;
		y = -y;
		z = -z;
	}

	/**
	* @return <code>true</code> if the given Vector3D is equivalent to this Vector3D.
	*/
	public boolean equals(Object a) {
		if (!(a instanceof Vector3D))
			return (false);
		return ((x == ((Vector3D) a).x) && (x == ((Vector3D) a).y) && (x == ((Vector3D) a).z));
	}

	public Vector3D mirror(Vector3D a) {
		return (this.add(a.mul(-2.0 * a.dot(this))));
	}

	public double getX() { return x; }
	public double getY() { return y; }
	public double getZ() { return z; }
	
	public double getCoordinate(int index) {
	
		switch (index) {
			case 0: return x;
			case 1: return y;
			case 2: return z;
			default: throw new IllegalArgumentException("illegal coordinates:" + index);
		}
	}

	/**
	*The modulus of this Vector3D
	* @return modulus
	*/
	public double modulus() {
		return (Math.sqrt(x * x + y * y + z * z));
	}

	/**
	* The modulus squared. <br>
	* @return modulus^2
	*/
	public double mod2() {
	return (x * x + y * y + z * z);
	}

	/**
	* Returns a new Vector3D that is this Vector3D multiplied by a scalar.
	*/
	public Vector3D mul(double a) {
	return (new Vector3D(x * a, y * a, z * a));
	}

	/**
	* Multiplies this Vector3D by a scalar.
	* @param a scalare da usare
	* @return questo stesso oggetto <code>Vector3D</code>
	*/
	public Vector3D mulU(double a) {
	x *= a;
	y *= a;
	z *= a;
	return (this);
	}

	/**
	* Subtracts the given Vector3D from this and returns the result.
	*/
	public Vector3D sub(Vector3D a) {
	return (new Vector3D(x - a.x, y - a.y, z - a.z));
	}

	/**
	* Subtracts the given Vector3D from this one.
	*/
	public Vector3D subU(Vector3D a) {
	x -= a.x;
	y -= a.y;
	z -= a.z;
	return (this);
	}

	/**
	* Returns this Vector3D with unit length.
	*/
	public Vector3D normalized() {
	return (this.mul(1.0 / this.modulus()));
	}

	/**
	* Transforms this Vector3D to be unit length.
	*/
	public Vector3D normalize() {
	double im = 1.0 / this.modulus();
	this.x *= im;
	this.y *= im;
	this.z *= im;
	return (this);
	}

	public String toString() {
	return ("Vector3D[" + x + "," + y + "," + z + "]");
	}

}