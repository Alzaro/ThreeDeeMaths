////////////////////////////////////////////////////////////////////////////////////////////////
// File :			EngineMath.cpp
// Author :			Corey Feist
// Purpose :		To define the EngineMath class which is a compilation of common 3D and 2D math functions
////////////////////////////////////////////////////////////////////////////////////////////////

#include "EngineMath.h"

//////////////////////////////////////////////////////////////////////////
// Common math functions
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
// General Utility functions
//////////////////////////////////////////////////////////////////////////

// Are two floating point numbers equal to each other
// Floating Point Error Safe
//
// IN:		a		The first number
//			b		The second number
//
// RETURN: TRUE iff |a-b| < Tolerance
//
// NOTE:	EPSILON is tolerance
bool IsEqual(float a, float b)
{
	// NOTE: Do not modify.
	return fabs(a - b) < EPSILON;
}

// Is a floating point value equal to zero
// Floating Point Error Safe
//
// IN:		a		The number to check
//
// RETURN:	TRUE iff |a| < Tolerance
//
// NOTE:	Tolerance set by EPSILON
bool IsZero(float a)
{
	// NOTE: Do not modify
	return (fabs(a))<EPSILON;
}

// RETURN: MAX of two numbers
float Max(float a, float b)
{
	// NOTE: Do not modify.
	return (a > b) ? a : b;
}

// RETURN: MIN of two numbers
float Min(float a, float b)
{
	// NOTE: Do not modify.
	return (a < b) ? a : b;
}

// RETURN: Converts input to radian measure
float Degrees_To_Radians(float Deg)
{
	// NOTE: Do not modify.
	return Deg * PI / 180.0f;
}

// RETURN: Converts input to degree measure
float Radians_To_Degrees(float Rad)
{
	// NOTE: Do not modify.
	return Rad * 180.0f / PI;
}
////////////////////////////////////////////////////////////////////////
// Linear Algebra Functions Day 1
///////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
// Vector Functions
//////////////////////////////////////////////////////////////////////////

// Check if two TVECTOR's are equal to each other
//
// IN:		v		First Vector
//			w		Second Vector
//
// RETURN:  True if v==w, False otherwise
//
// NOTE:	Use's all four components
//			Should be floating point error safe.
bool Vector_IsEqual(TVECTOR v, TVECTOR w)
{
	if (IsEqual(v.x,w.x) && IsEqual(v.y,w.y) && IsEqual(v.z,w.z) && IsEqual(v.w,w.w))
	{
		return true;
	}
	return false;
}

// ADD two TVECTOR's togother
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v + w
//
// NOTE:	Use's all four components
TVECTOR Vector_Add(TVECTOR v, TVECTOR w)
{
	v.x = v.x + w.x;
	v.y = v.y + w.y;
	v.z = v.z + w.z;
	v.w = v.w + w.w;
	return v;
}

// SUBTRACT one TVECTOR from another
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v - w
//
// NOTE:	Use's all four components
TVECTOR Vector_Sub(TVECTOR v, TVECTOR w)
{
	v.x = v.x - w.x;
	v.y = v.y - w.y;
	v.z = v.z - w.z;
	v.w = v.w - w.w;
	return v;
}

// MULTIPLY all four components of a TVECTOR by a scalar
//
// IN:		v		The vector to scale
//			s		The value to scale by
//
// RETURN:  s * v
TVECTOR Vector_Scalar_Multiply(TVECTOR v, float s)
{
	v.x *= s;
	v.y *= s;
	v.z *= s;
	v.w *= s;
	return v;
}

// NEGATE all the components of a TVECTOR
//
// IN:		v		The vector to negate
//
// RETURN:	-1 * v
//
// NOTE:	Use's all four components
TVECTOR Vector_Negate(TVECTOR v)
{

	TVECTOR temp;
	for (int i = 0; i < 4; ++i)
	{
		temp.e[i] = v.e[i] * -1;
	}
	return temp;
}

// Perform a Dot Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (DOT) w
//
// NOTE:	Use's all four components
float Vector_Dot(TVECTOR v, TVECTOR w)
{
	float result = 0;
	for (int i = 0; i < 4; ++i)
	{
		result += (v.e[i] * w.e[i]);
	}
	return result;
}

// Perform a Cross Product on two TVECTOR's
//
// IN:		v		First Vector. Left Hand Side
//			w		Second Vector. Right Hand Side
//
// RETURN:  v (CROSS) w
//
// NOTE:	The w-component of each vector is not used.
//			The resultant vector will have a w-component of zero.
TVECTOR Vector_Cross(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	temp.x = (v.y * w.z) - (v.z * w.y);
	temp.y = -((v.x * w.z) - (v.z * w.x));
	temp.z = (v.x * w.y) - (v.y * w.x);
	temp.w = 0;
	return temp;
}

// Find the squared length of a TVECTOR
//
// IN:		v		The vector to find the squared length of
//
// RETURN:	Squared Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_LengthSq(TVECTOR v)
{
	float result = 0;
	for (int i = 0; i < 4; ++i)
	{
		result += v.e[i] * v.e[i];
	}
	return result;
}

// Find the length of a TVECTOR
//
// IN:		v		The vector to find the length of
//
// RETURN:	Length of TVECTOR
//
// NOTE:	Use's all four components
float Vector_Length(TVECTOR v)
{
	float result = 0;
	for (int i = 0; i < 4; ++i)
	{
		result += v.e[i] * v.e[i];
	}
	result = sqrt(result);
	return result;
}

// Normalize a TVECTOR
//
// IN:		v		The vector to normalize
//
// RETURN:	Normalized version of v
//
// NOTE:	Use's all four components
TVECTOR Vector_Normalize(TVECTOR v)
{
	TVECTOR temp;
	for (int i = 0; i < 4; ++i)
	{
		temp.e[i] = (v.e[i] / Vector_Length(v));
	}
	return temp;
}

// Makes a TVECTOR's w-component normalized
//
// IN:		v		The vector (point object) to homogenise
//
// RETURN:	The homogenised vector (point)
//
// NOTE:	If the w-component of the vector is 0 then the
//			function will return a zero vector with a w-component
//			of 0.
TVECTOR Vector_Homogenise(TVECTOR v)
{
	TVECTOR temp;
	if (v.w == 0)
	{
		for (int i = 0; i < 4; ++i)
		{
			temp.e[i] = 0;
		}
		return temp;
	}
	else
	{
		for (int i = 0; i < 4; ++i)
		{
			temp.e[i] =( v.e[i] / v.w);
		}
	}
	return temp;
}

// Get a TVECTOR made from the maximun components of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A maximized vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Maximize(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	for (int i = 0; i < 4; ++i)
	{
		if (v.e[i] > w.e[i])
		{
			temp.e[i] = v.e[i];
		}
		else 
		{
			temp.e[i] = w.e[i];
		}
	}
	return temp;
}

// Get a TVECTOR made from the minimum components of two TVECTOR's
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A minimum vector
//
// NOTE:	Use's all four components
TVECTOR Vector_Minimize(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	for (int i = 0; i < 4; ++i)
	{
		if (v.e[i] < w.e[i])
		{
			temp.e[i] = v.e[i];
		}
		else
		{
			temp.e[i] = w.e[i];
		}
	}
	return temp;
}

// Get a TVECTOR made from the average of two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:	A vector made from the average of two vectors
//
// NOTE:	Use's all four components

TVECTOR Vector_Average(TVECTOR v, TVECTOR w)
{
	TVECTOR temp;
	int i = 0;
	for (; i < 4; ++i)
	{
		temp.e[i] = v.e[i] + w.e[i];
	}
	
	for (i = 0; i < 4; ++i)
	{
		temp.e[i] = temp.e[i] * 0.5;
	}
	return temp;
}

// Find the angle between two TVECTORs
//
// IN:		v		The first vector
//			w		The second vector
//
// RETURN:  The angle in degrees between the two vectors
//
// NOTE:	If either vector is a zero vector then the return
//			value will be 0.
float Vector_AngleBetween(TVECTOR v, TVECTOR w)
{
	float result;
	result = acos((Vector_Dot(v, w) / (Vector_Length(v) * Vector_Length(w))));
	result = Radians_To_Degrees(result);
	return result;
}

// Get the distance one TVECTOR points in the direction of another
// TVECTOR
//
// IN:		v		The first vector
//			w		The direction of the component
//
// RETURN:	The distance that v points in the direction of w.
//
// NOTE:	If w or v is a zero vector then the return value is zero.
float Vector_Component(TVECTOR v, TVECTOR w)
{
	float result;
	result = (Vector_Dot(v, w) / Vector_Length(w));
	return result;
}

// Get the TVECTOR that represents v projected on w.
//
// IN:		v		The first vector
//			w		The direction of the projection
//
// RETURN:	The projection of v onto w
//
// NOTE:	If w or v is a zero vector then the return value is zero.
TVECTOR Vector_Project(TVECTOR v, TVECTOR w)
{
	if ((v.x == 0 && v.y == 0 && v.z == 0 && v.w == 0) || (w.x == 0 && w.y == 0 && w.z == 0 && w.w == 0))
	{
		TVECTOR temp;
		for (int i = 0; i < 4; ++i)
		{
			temp.e[i] = 0;
		}
		return temp;
	}
	TVECTOR temp;
	TVECTOR normVec = Vector_Normalize(w);
	float comp = Vector_Component(v, w);
	for (int i = 0; i < 4; ++i)
	{
		temp.e[i] = comp * normVec.e[i];
	}
	return temp;
}


// Get the reflection of v across w
//
// IN:		v		The vector to reflect
//			w		The "axis" to reflect across
//
// RETURN:	v reflected across w
//
// NOTE:	If w is a zero vector then return -v.
TVECTOR Vector_Reflect(TVECTOR v, TVECTOR w)
{
	if (w.x == 0 && w.y == 0 && w.z == 0 && w.w == 0)
	{
		return Vector_Negate(v);
	}
	TVECTOR temp;
	TVECTOR norm = Vector_Normalize(w);
	float dot = Vector_Dot(v, norm);
	for (int i = 0; i < 4; ++i)
	{
		temp.e[i] = v.e[i] - (2 * norm.e[i] * dot);
	}
	temp = Vector_Negate(temp);
	return temp;
}

//////////////////////////////////////////////////////////////////////////
// Matrix Functions
//////////////////////////////////////////////////////////////////////////

// Get a [0] matrix
//
// RETURN: A 0 4x4 matrix
TMATRIX Matrix_Zero(void)
{
	TMATRIX m;
	for (int i = 0; i < 16; ++i)
	{
		m.e[i] = 0;
	}

	return m;
}

// Get a [I] matrix
//
// RETURN: A 4x4 Identity matrix
TMATRIX Matrix_Identity(void)
{
	TMATRIX m{ 1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1 };
	return m;
}

// Get a translation matrix
//
// IN:		x		Amount of translation in the x direction
//			y		Amount of translation in the y direction
//			z		Amount of translation in the z direction
//
// RETURN:	The translation matrix
TMATRIX Matrix_Create_Translation(float x, float y, float z)
{
	TMATRIX m = { 1,0,0,x,
				  0,1,0,y,
				  0,0,1,z,
				  0,0,0,1 };
	return m;
}

// Create a scale matrix
//
// IN:		x		Amount to scale in the x direction
//			y		Amount to scale in the y direction
//			z		Amount to scale in the z direction
//
// RETURN:	The scale matrix
TMATRIX Matrix_Create_Scale(float x, float y, float z)
{
	TMATRIX m = { x,0,0,0,
			      0,y,0,0,
			      0,0,z,0,
			      0,0,0,1 };
	return m;
}

// Get a rotation matrix for rotation about the x-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A X-Rotation Matrix
TMATRIX Matrix_Create_Rotation_X(float Deg)
{
	TMATRIX m = { 1,0,0,0,
				  0,cos(Degrees_To_Radians(Deg)),sin(Degrees_To_Radians(Deg)) * -1,0,
				  0,sin(Degrees_To_Radians(Deg)),cos(Degrees_To_Radians(Deg)),0,
				  0,0,0,1 };
	return m;
}

// Get a rotation matrix for rotation about the y-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Y-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Y(float Deg)
{
	TMATRIX m = { cos(Degrees_To_Radians(Deg)),0,sin(Degrees_To_Radians(Deg)),0,
				  0,1,0,0,
				  sin(Degrees_To_Radians(Deg)) * -1,0,cos(Degrees_To_Radians(Deg)),0,
				  0,0,0,1 };
	return m;
}

// Get a rotation matrix for rotation about the z-axis
//
// IN:		Deg		Angle to rotate ( Degree measure)
//
// RETURN:	A Z-Rotation Matrix
TMATRIX Matrix_Create_Rotation_Z(float Deg)
{
	TMATRIX m = { cos(Degrees_To_Radians(Deg)),sin(Degrees_To_Radians(Deg)) * -1,0,0,
			      sin(Degrees_To_Radians(Deg)),cos(Degrees_To_Radians(Deg)),0,0,
				  0,0,1,0,
				  0,0,0,1 };
	return m;
}

// ADD two matrices together
//
// IN:		m		The first matrix
//			n		The second matrix
//
// RETURN: m + n
TMATRIX Matrix_Matrix_Add(TMATRIX m, TMATRIX n)
{
	TMATRIX temp;
	for (int i = 0; i < 16; ++i)
	{
		temp.e[i] = m.e[i] + n.e[i];
	}
	return temp;
}

// SUBTRACT two matrices
//
// IN:		m		The first matrix (left hand side)
//			n		The second matrix (right hand side)
//
// RETURN: m - n
TMATRIX Matrix_Matrix_Sub(TMATRIX m, TMATRIX n)
{
	TMATRIX temp;
	for (int i = 0; i < 16; ++i)
	{
		temp.e[i] = m.e[i] - n.e[i];
	}
	return temp;
}

// Multiply a matrix by a scalar
//
// IN:		m		The matrix to be scaled (right hand side)
//			s		The value to scale by   (left hand side)
//
// RETURN:	The matrix formed by s*[m]
TMATRIX Matrix_Scalar_Multiply(TMATRIX m, float s)
{
	TMATRIX temp;
	for (int i = 0; i < 16; ++i)
	{
		temp.e[i] = m.e[i] * s;
	}
	return temp;
}

// Negate a matrix
//
// IN:		m		The matrix to negate
//
// RETURN:  The negation of m
TMATRIX Matrix_Negate(TMATRIX m)
{
	TMATRIX temp;
	for (int i = 0; i < 16; ++i)
	{
		temp.e[i] = m.e[i] * -1;
	}
	return temp;
}

// Transpose a matrix
//
// IN:		m		The matrix to transpose
//
// RETURN:	The transpose of m
TMATRIX Matrix_Transpose(TMATRIX m)
{
	TMATRIX temp;
	temp._e11 = m._e11;
	temp._e12 = m._e21;
	temp._e13 = m._e31;
	temp._e14 = m._e41;
	temp._e21 = m._e12;
	temp._e22 = m._e22;
	temp._e23 = m._e32;
	temp._e24 = m._e42;
	temp._e31 = m._e13;
	temp._e32 = m._e23;
	temp._e33 = m._e33;
	temp._e34 = m._e43;
	temp._e41 = m._e14;
	temp._e42 = m._e24;
	temp._e43 = m._e34;
	temp._e44 = m._e44;
	return temp;
}

// Multipy a matrix and a vector
//
// IN:		m		The matrix (left hand side)
//			v		The vector (right hand side)
//
// RETURN:	[m]*v
TVECTOR Matrix_Vector_Multiply(TMATRIX m, TVECTOR v)
{
	TVECTOR temp;
	temp.x = m._e11 * v.x + m._e12 * v.y + m._e13 * v.z + m._e14 * v.w;
	temp.y = m._e21 * v.x + m._e22 * v.y + m._e23 * v.z + m._e24 * v.w;
	temp.z = m._e31 * v.x + m._e32 * v.y + m._e33 * v.z + m._e34 * v.w;
	temp.w = m._e41 * v.x + m._e42 * v.y + m._e43 * v.z + m._e44 * v.w;
	return temp;
}

// Multipy a vector and a matrix
//
// IN:		v		The vector ( left hand side)
//			m		The matrix (right hand side)
//
// RETURN:	v*[m]
TVECTOR Vector_Matrix_Multiply(TVECTOR v, TMATRIX m)
{
	TVECTOR temp;
	temp.x = m._e11 * v.x + m._e21 * v.y + m._e31 * v.z + m._e41 * v.w;
	temp.y = m._e12 * v.x + m._e22 * v.y + m._e32 * v.z + m._e42 * v.w;
	temp.z = m._e13 * v.x + m._e23 * v.y + m._e33 * v.z + m._e43 * v.w;
	temp.w = m._e14 * v.x + m._e24 * v.y + m._e34 * v.z + m._e44 * v.w;
	return temp;
}
// Multiply a matrix by a matrix
//
// IN:		m		First Matrix (left hand side)
//			n		Second Matrix (right hand side)
//
// RETURN:	[m]*[n]
TMATRIX Matrix_Matrix_Multiply(TMATRIX m, TMATRIX n)
{
	TMATRIX temp;						 
	temp._e11 = m._e11 * n._e11 + m._e12 * n._e21 + m._e13 * n._e31 + m._e14 * n._e41;
	temp._e12 = m._e11 * n._e12 + m._e12 * n._e22 + m._e13 * n._e32 + m._e14 * n._e42;
	temp._e13 = m._e11 * n._e13 + m._e12 * n._e23 + m._e13 * n._e33 + m._e14 * n._e43;
	temp._e14 = m._e11 * n._e14 + m._e12 * n._e24 + m._e13 * n._e34 + m._e14 * n._e44;
	temp._e21 = m._e21 * n._e11 + m._e22 * n._e21 + m._e23 * n._e31 + m._e24 * n._e41;
	temp._e22 = m._e21 * n._e12 + m._e22 * n._e22 + m._e23 * n._e32 + m._e24 * n._e42;
	temp._e23 = m._e21 * n._e13 + m._e22 * n._e23 + m._e23 * n._e33 + m._e24 * n._e43;
	temp._e24 = m._e21 * n._e14 + m._e22 * n._e24 + m._e23 * n._e34 + m._e24 * n._e44;
	temp._e31 = m._e31 * n._e11 + m._e32 * n._e21 + m._e33 * n._e31 + m._e34 * n._e41;
	temp._e32 = m._e31 * n._e12 + m._e32 * n._e22 + m._e33 * n._e32 + m._e34 * n._e42;
	temp._e33 = m._e31 * n._e13 + m._e32 * n._e23 + m._e33 * n._e33 + m._e34 * n._e43;
	temp._e34 = m._e31 * n._e14 + m._e32 * n._e24 + m._e33 * n._e34 + m._e34 * n._e44;
	temp._e41 = m._e41 * n._e11 + m._e42 * n._e21 + m._e43 * n._e31 + m._e44 * n._e41;
	temp._e42 = m._e41 * n._e12 + m._e42 * n._e22 + m._e43 * n._e32 + m._e44 * n._e42;
	temp._e43 = m._e41 * n._e13 + m._e42 * n._e23 + m._e43 * n._e33 + m._e44 * n._e43;
	temp._e44 = m._e41 * n._e14 + m._e42 * n._e24 + m._e43 * n._e34 + m._e44 * n._e44;
	return temp;
}

////////////////////////////////////////////////////////////////////////
// Matrix Functions
///////////////////////////////////////////////////////////////////////

// HELPER FUNCTION  *** NOT GRADED, ONLY SUGGESTED ***
// USE THIS FUNCTION TO FIND THE DETERMINANT OF A 3*3
// MATRIX. IT CAN BE USED IN THE MATRIX DETERMINANT
// AND MATRIX INVERSE FUNCTIONS BELOW
// 
// RETURN:	The determinant of a 3x3 matrix
float Matrix_DeterminantHELP(float e_11,float e_12,float e_13,
						     float e_21,float e_22,float e_23,
						     float e_31,float e_32,float e_33)
{
	float result, a, b, c;
	a = ((e_22 * e_33) - (e_23 * e_32)) * e_11;
	b = ((e_21 * e_33) - (e_23 * e_31)) * e_12;
	c = ((e_21 * e_32) - (e_22 * e_31)) * e_13;
	result = a + (b*-1) + c;
	return result;
}

// Get the determinant of a matrix
//
// IN:		m		The ONE!
//
// RETURN:	It's deterinant
float Matrix_Determinant(TMATRIX m)
{
	float det = 0;
	det += m._e11 * Matrix_DeterminantHELP(m._e22, m._e23, m._e24, m._e32, m._e33, m._e34, m._e42, m._e43, m._e44);//+
	det -= m._e12 * Matrix_DeterminantHELP(m._e21, m._e23, m._e24, m._e31, m._e33, m._e34, m._e41, m._e43, m._e44);//-
	det += m._e13 * Matrix_DeterminantHELP(m._e21, m._e22, m._e24, m._e31, m._e32, m._e34, m._e41, m._e42, m._e44);//+
	det -= m._e14 * Matrix_DeterminantHELP(m._e21, m._e22, m._e23, m._e31, m._e32, m._e33, m._e41, m._e42, m._e43);//-
	return det;
}

// Get the inverse of a matrix
//
// IN:		m		The matrix to inverse
//
// RETURN:	The Inverse of [m]
//
// NOTE: Returns the matrix itself if m is not invertable.
TMATRIX Matrix_Inverse(TMATRIX m)
{
	
	if (Matrix_Determinant(m) != 0)
	{
		TMATRIX tran;
		TMATRIX cofactor{ m._e11, -m._e12, m._e13, -m._e14,
						 -m._e21, m._e22, -m._e23, m._e24,
						 m._e31, -m._e32, m._e33, -m._e34,
						 -m._e41, m._e42, -m._e43, m._e44 };
		tran = Matrix_Transpose(m);
		float detCoFactor = Matrix_Determinant(cofactor);
		float det = 1 / detCoFactor;
		for (int i = 0; i < 16; ++i)
		{
			cofactor.e[i] *= det;
		}
		return cofactor;
	}
	return m;
}

