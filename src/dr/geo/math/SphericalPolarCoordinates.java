/*
 * SphericalPolarCoordinates.java
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

import dr.evolution.continuous.Contrastable;

public class SphericalPolarCoordinates implements Contrastable {
	
	/**
	 * Create spherical polar coordinates from the given vector
	 */
	public SphericalPolarCoordinates(Vector3D v) {
		this(v, v.modulus());
	}
	
	/**
	 * Create spherical polar coordinates from the given vector
	 */
	public SphericalPolarCoordinates(Vector3D v, double newRadius) {
	
		radius = v.modulus();
		theta = Math.acos(v.getZ()/radius);
		phi = Math.atan(v.getY()/v.getX());
		if (v.getY() < 0) phi += Math.PI;
		
		radius = newRadius;
	}
	
	/**
	 * Create spherical polar coordinates from given latitude and longitude
	 */
	public SphericalPolarCoordinates(double latitude, double longitude) {
        this(latitude,longitude,VOLUMETRIC_RADIUS_OF_EARTH);
    }

    public SphericalPolarCoordinates(double latitude, double longitude, double radius) {
	
		this.radius = radius;
		theta = (90.0-latitude)*Math.PI/180.0;
		
		if (longitude < 0) longitude += 360;
		
		phi = longitude * Math.PI/180;
	}
	
	public Vector3D getCartesianCoordinates() {
		double x = radius * Math.sin(theta) * Math.cos(phi);
		double y = radius * Math.sin(theta) * Math.sin(phi);
		double z = radius * Math.cos(theta);
	
		return new Vector3D(x, y, z);
	}
	
	/**
	 * @return the latitude of these spherical polar coordinates in degrees.
	 */
	public double getLatitude() {
		return 90 - (theta * 180.0 / Math.PI);
	}
	
	/**
	 * @return the longitude of these spherical polar coordinates in degrees.
	 */
	public double getLongitude() {
		double degrees = phi * 180.0 / Math.PI;
		if (degrees > 180.0) {
			degrees -= 360;
		}
		return degrees;
	}
	
	/**
	 * @return the angle (along the great circle) between this position and the given position.
	 */
	public double angle(SphericalPolarCoordinates B) {
		
		Vector3D v1 = getCartesianCoordinates().normalized();
		Vector3D v2 = B.getCartesianCoordinates().normalized();
		
		double sinAngle = v1.cross(v2).modulus();
		
		double cosAngle = v1.dot(v2);
		
		double angle = Math.atan2(sinAngle, cosAngle);
		
		return angle;
	}
	
	/**
	 * @return the great circle distance between this position and the given position.
	 */
	public double distance(SphericalPolarCoordinates B) {
		
		return angle(B) * radius;
	}
	
	/**
	 * @return the mid point of the two positions.
	 */
	public static SphericalPolarCoordinates getMidpoint(SphericalPolarCoordinates A, SphericalPolarCoordinates B) {
		
		return interpolate(A, B, 0.5);
	
	}
	
	/**
	 * alpha = 0 implies returning A, alpha = 1 implies returning B
	 * @return a point between the two positions.
	 */
	public static SphericalPolarCoordinates interpolate(SphericalPolarCoordinates A, SphericalPolarCoordinates B, double alpha) {
		
		Vector3D from = A.getCartesianCoordinates().normalized();
		Vector3D to = B.getCartesianCoordinates().normalized();
		
		double cosang = from.dot(to);
  		if (cosang < 0) {
  			cosang = -cosang;
  			to.negate();
  			System.out.println("Argg");
  		}

		double fracFrom, fracTo;

		if (cosang > (0.999995) ) {  // small angle limit
			fracFrom = alpha;
			fracTo   = 1.0-alpha;
		} else {  // normal case
			double ang = Math.acos( cosang );
			fracFrom = Math.sin( alpha * ang);
			fracTo   = Math.sin( (1.0-alpha) * ang );
		}
		
		Vector3D inbetween = from.mul(fracFrom).add(to.mul(fracTo)).normalized();
		
		return new SphericalPolarCoordinates(inbetween, A.getRadius());
	}
	
	public double getRadius() { return radius; }
	
	public String toString() {

        return toLatLongString();
        //return "polarCoordinates[theta="+theta+", phi=" + phi + ", r=" + radius + "]";
	}
	
	public String toLatLongString() {
	
		double lat1 = getLatitude();
		double long1 = getLongitude();
		
		String latSuffix = " N";
		if (lat1 < 0) {
			lat1 = -lat1;
			latSuffix = " S";
		}
		
		String longSuffix = " E";
		if (long1 < 0) {
			long1 = -long1;
			longSuffix = " W";
		}
		
		return lat1 + latSuffix + ", " + long1 + longSuffix;
	}
	
	public static void main(String[] args) {
	
		double latitude1 = 66.0;
		double longitude1 = -174.0;
	
		double latitude2 = 64.0;
		double longitude2 = -142.0;
	
		SphericalPolarCoordinates p1 = new SphericalPolarCoordinates(latitude1, longitude1);
		SphericalPolarCoordinates p2 = new SphericalPolarCoordinates(latitude2, longitude2);
		
		SphericalPolarCoordinates p1a = new SphericalPolarCoordinates(p1.getCartesianCoordinates());
		SphericalPolarCoordinates p2a = new SphericalPolarCoordinates(p2.getCartesianCoordinates());
		
		
		System.out.println("point1="+p1.toLatLongString());
		System.out.println("point2="+p2.toLatLongString());
		
		System.out.println("point1_test="+p1a.toLatLongString());
		System.out.println("point2_test="+p2a.toLatLongString());
		
		
		System.out.println("distance=" + p1.distance(p2) + " kilometres");
		
		System.out.println("midpoint=" + getMidpoint(p1, p2).toLatLongString());
	}



    //************************************************************************
	// Contrastable interface
	//************************************************************************
	
	public double getDifference(Contrastable spc) {
	
		if (spc instanceof SphericalPolarCoordinates) {
			return distance((SphericalPolarCoordinates)spc);
		} else throw new IllegalArgumentException("Expected a spherical polar coordinate");
	}
	
	public Contrastable getWeightedMean(double weight1, Contrastable spc1, double weight2, Contrastable spc2) {
		
		double alpha = weight1 / (weight1 + weight2);
		
		if (spc1 instanceof SphericalPolarCoordinates && spc2 instanceof SphericalPolarCoordinates) {
			return interpolate((SphericalPolarCoordinates)spc1, (SphericalPolarCoordinates)spc2, alpha);
		} else throw new IllegalArgumentException("Expected a spherical polar coordinate");
	}

	//************************************************************************
	// Private members
	//************************************************************************
	
	// angle off of north pole (0 to pi)
	private double theta = 0.0;
	
	// angle east of greenwich (0 to 2*pi)
	private double phi = 0.0;	
	
	// radius
	private double radius = 1.0;
	
	private static double VOLUMETRIC_RADIUS_OF_EARTH = 6371.0;
}