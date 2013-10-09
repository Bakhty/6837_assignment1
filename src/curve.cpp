#include "curve.h"
#include "extra.h"
using namespace std;

namespace
{
    // Approximately equal to.  We don't want to use == because of
    // precision issues with floating point.
    inline bool approx( const Vector3f& lhs, const Vector3f& rhs )
    {
        const float eps = 1e-8f;
        return ( lhs - rhs ).absSquared() < eps;
    }

    
}

const Matrix4f Bernstein(1.0f,-3.0f,3.0f,-1.0f,
    0.0f,3.0f,-6.0f,3.0f,
    0.0f,0.0f,3.0f,-3.0f,
    0.0f,0.0f,0.0f,1.0f);
const Matrix4f dBernstein(-3,6,-3,0,3,-12,9,0,0,6,-9,0,0,0,3,0);
const Matrix4f BernsteinInv(1,1,1,1,0,1.f/3,2.f/3,1,0,0,1.f/3,1,0,0,0,1);
const Matrix4f BSpline(1.f/6,-1.f/2,1.f/2,-1.f/6,
    2.f/3,0,-1,1.f/2,
    1.f/6,1.f/2,1.f/2,-1.f/2,
    0,0,0,1.f/6);

Curve evalBezier( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 || P.size() % 3 != 1 )
    {
        cerr << "evalBezier must be called with 3n+1 control points." << endl;
        exit( 0 );
    }

    // You should implement this function so that it returns a Curve
    // (e.g., a vector< CurvePoint >).  The variable "steps" tells you
    // the number of points to generate on each piece of the spline.
    // At least, that's how the sample solution is implemented and how
    // the SWP files are written.  But you are free to interpret this
    // variable however you want, so long as you can control the
    // "resolution" of the discretized spline curve with it.

    // Also note that you may assume that all Bezier curves that you
    // receive have G1 continuity.  Otherwise, the TNB will not be
    // be defined at points where this does not hold.

    cerr << "\t>>> evalBezier has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    // initialize curve
    Curve C;

    // iterate every four control points
    for (unsigned i = 0; i < P.size() - 1; i += 3) {
        // Geometry matrix
        Matrix4f G(Vector4f(P[i],0), Vector4f(P[i+1],0), Vector4f(P[i+2],0), Vector4f(P[i+3],0));
        
        // for use in inner for loop
        // Geometry * Spline Basis constants
        Matrix4f GxB = G*Bernstein;
        Matrix4f GxdB = G*dBernstein;

        for (unsigned s = 0; s <= steps; s++) {
            // time step
            float t = s*1.0/steps;
            // Power basis
            Vector4f basis(1,t,pow(t,2),pow(t,3));

            // V = Q(t) = GBT(t)
            Vector3f V = (GxB*basis).xyz();
            // T = Q'(t).normalized()
            Vector3f T = (GxdB*basis).xyz().normalized();
            Vector3f N;
            Vector3f B;

            // Ni = (Bi-1xTi).normalized()
            // Bi = (TixNi).normalized()
            if (i == 0 and s == 0) {
                B = Vector3f(0,0,1); // put normals on left
                // make sure T and B are not parallel
                if (Vector3f::cross(T,B) == Vector3f(0,0,0)) {
                    B = Vector3f(0,1,0);
                }
                N = Vector3f::cross(B,T).normalized();
            } else {
                CurvePoint last = C.back();
                N = Vector3f::cross(last.B,T).normalized();
                B = Vector3f::cross(T,N).normalized();
            }

            CurvePoint next = {V,T,N,B};
            C.push_back(next);

        }

    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning Bezier curve." << endl;

    return C;
}

Curve evalBspline( const vector< Vector3f >& P, unsigned steps )
{
    // Check
    if( P.size() < 4 )
    {
        cerr << "evalBspline must be called with 4 or more control points." << endl;
        exit( 0 );
    }

    cerr << "\t>>> evalBSpline has been called with the following input:" << endl;

    cerr << "\t>>> Control points (type vector< Vector3f >): "<< endl;
    for( unsigned i = 0; i < P.size(); ++i )
    {
        cerr << "\t>>> " << P[i] << endl;
    }

    // Change basis from B-spline to Bezier then call evalBezier

    vector<Vector3f> bsplineP;
    for (unsigned i = 0; i < P.size()-3; i++) {
        // Geometry matrix, add row of 0s to multiply by other 4x4s
        Matrix4f G(Vector4f(P[i],0), Vector4f(P[i+1],0), Vector4f(P[i+2],0), Vector4f(P[i+3],0));
        // make new control points G*B1*B2_inv
        Matrix4f newG = G*BSpline*BernsteinInv;
        // iterate through columns of new basis G to get the new vectors
        for (unsigned j = (i==0 ? 0 : 1); j < 4; j++) {
            Vector3f newPoint = newG.getCol(j).xyz();
            bsplineP.push_back(newPoint);
        }
    }

    cerr << "\t>>> Steps (type steps): " << steps << endl;
    cerr << "\t>>> Returning B-spline curve." << endl;

    Curve bsplinecurve = evalBezier(bsplineP, steps);

    return bsplinecurve;
}

Curve evalCircle( float radius, unsigned steps )
{
    // This is a sample function on how to properly initialize a Curve
    // (which is a vector< CurvePoint >).
    
    // Preallocate a curve with steps+1 CurvePoints
    Curve R( steps+1 );

    // Fill it in counterclockwise
    for( unsigned i = 0; i <= steps; ++i )
    {
        // step from 0 to 2pi
        float t = 2.0f * M_PI * float( i ) / steps;

        // Initialize position
        // We're pivoting counterclockwise around the y-axis
        R[i].V = radius * Vector3f( cos(t), sin(t), 0 );
        
        // Tangent vector is first derivative
        R[i].T = Vector3f( -sin(t), cos(t), 0 );
        
        // Normal vector is second derivative
        R[i].N = Vector3f( -cos(t), -sin(t), 0 );

        // Finally, binormal is facing up.
        R[i].B = Vector3f( 0, 0, 1 );
    }

    return R;
}

void drawCurve( const Curve& curve, float framesize )
{
    // Save current state of OpenGL
    glPushAttrib( GL_ALL_ATTRIB_BITS );

    // Setup for line drawing
    glDisable( GL_LIGHTING ); 
    glColor4f( 1, 1, 1, 1 );
    glLineWidth( 1 );
    
    // Draw curve
    glBegin( GL_LINE_STRIP );
    for( unsigned i = 0; i < curve.size(); ++i )
    {
        glVertex( curve[ i ].V );
    }
    glEnd();

    glLineWidth( 1 );

    // Draw coordinate frames if framesize nonzero
    if( framesize != 0.0f )
    {
        Matrix4f M;

        for( unsigned i = 0; i < curve.size(); ++i )
        {
            M.setCol( 0, Vector4f( curve[i].N, 0 ) );
            M.setCol( 1, Vector4f( curve[i].B, 0 ) );
            M.setCol( 2, Vector4f( curve[i].T, 0 ) );
            M.setCol( 3, Vector4f( curve[i].V, 1 ) );

            glPushMatrix();
            glMultMatrixf( M );
            glScaled( framesize, framesize, framesize );
            glBegin( GL_LINES );
            glColor3f( 1, 0, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 1, 0, 0 );
            glColor3f( 0, 1, 0 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 1, 0 );
            glColor3f( 0, 0, 1 ); glVertex3d( 0, 0, 0 ); glVertex3d( 0, 0, 1 );
            glEnd();
            glPopMatrix();
        }
    }
    
    // Pop state
    glPopAttrib();
}

