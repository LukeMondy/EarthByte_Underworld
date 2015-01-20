/*
	 File contains utility functions for Underworld Spherical

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <float.h>

#include <petsc.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>
#include <Underworld/Underworld.h>

#include "types.h"
#include "Components.h"
/*
   How the mesh is specified

         Y
         |
         |
         |  /
         | /
         |/
         |----------- X
        /
       /
      /
     Z

   The domain is specified in spherical coordinates, (r,longitude, latitude).
   To define these spaces we use the existing nomencluture:

   elementResI = discretisation of radius, the radial length  range (0, inf)
   elementResJ = discretisation of longitude, angle in the X-Z plane - 0 along X-axis. 
                 RH screw rule thumb pointing along Y-axis (-180,180)
   elementResK = discretisation of latitude, the angle in the Y-Z plane - 0 is along Z-axis. (-90, 90)


   Implementation notes:
   All spherical geometry calculation are performed in radians to match native c functions. 
   As input geometry is in in degrees so should output geometric data, to facilitate Spherical_Radians2Output()

*/

void Spherical_Radians2Output( double *input ) {
   /* 
      Converts input angles, in RADIANS, to output angles in DEGREES.
      Assumed vectors is (radius, longitude, latitude)
   */
   input[1] = 180/M_PI * input[1];
   input[2] = 180/M_PI * input[2];
}

void Spherical_GetRotationMatrixIJK_RSNodes( Mesh* mesh, int dNodeID, double* rot )
{

   /* Returns the rotation matrix for a boundary node for the RegionalSphere mesh (RSGenerator)

      Canonoical mesh coords r, eta, xi
      -45 <= eta, xi <= 45

      s = (1+ tan(eta)*tan(eta) + tan(xi)*tan(xi) )^(-0.5)
      
      x = r * s * tan(eta) 
      y = r * s * tan(xi) 
      z = r * s

      trick is to calculate the unit vectors for:
      the normal - n
      the tangent vector north - tn
      the tangent vector east - te

      then apply those vectors in the matrix below, which is the rotation matrix

      | dx_dt |     | n_x  tn_x  te_x |  | dr_dt   |
      | dy_dt |  =  | n_y  tn_y  te_y |  | deta_dt |
      | dz_dt |     | n_z  tn_z  te_z |  | dxi_dt  |

      To calculate the n, tn, te vectors we use the following definition.

      | dx_dt |     | dx_dr  dx_deta  dx_dxi |  | dr_dt   |
      | dy_dt |  =  | dy_dr  dy_deta  dy_dxi |  | deta_dt |
      | dz_dt |     | dz_dr  dz_deta  dz_dxi |  | dxi_dt  |

      1) Evaulate the coeffient matrix.
      2) each column in the matrix is either n, tn, te
         we convert these to unit vector. And this is the rotation matrix
   */
         
   ProjectionGenerator* gen = NULL;
   unsigned inds[3];
   double *vert, centre[3], vec[3], dimRes[3], R, x_z_len;
   double mag;
   double r, t, p;
   unsigned gNode;

   Grid **grid = (Grid** )Mesh_GetExtension(mesh,Grid**,mesh->vertGridId);
   centre[0] = centre[1] = centre[2] = 0.0;
   vert = Mesh_GetVertex( mesh, dNodeID );
   Grid_Lift( grid[0], dNodeID, inds );

   // calculate the radius of projection
   r = sqrt( vert[0]*vert[0] + vert[1]*vert[1] + vert[2]*vert[2] );
   assert( r>1e-9 ); // sanity check for zero radius

   double r_vec[3], tang1_vec[3], tang2_vec[3];
   double xi  = atan2(vert[0], vert[2] );
   double eta = atan2(vert[1], vert[2] );

   r_vec[0] = vert[0]/r;
   r_vec[1] = vert[1]/r;
   r_vec[2] = vert[2]/r;

   /*
   if( inds[0] == grid[0]->sizes[0]-1 ) {
      double t = xi;
      double p = asin( vert[1]/r );
      rot[0] = sin(t)*cos(p) ;  rot[1] = r*cos(t)*cos(p)  ;  rot[2] = -r*sin(t)*sin(p);
      rot[3] = sin(p)        ;  rot[4] = 0                ;  rot[5] = r*cos(p);
      rot[6] = cos(t)*cos(p) ;  rot[7] = -r*sin(t)*cos(p) ;  rot[8] = -r*cos(t)*sin(p);

      mag = sqrt( rot[0]*rot[0] + rot[3]*rot[3] + rot[6]*rot[6] );
      rot[0] = rot[0] / mag;
      rot[3] = rot[3] / mag;
      rot[6] = rot[6] / mag;

      mag = sqrt( rot[1]*rot[1] + rot[4]*rot[4] + rot[7]*rot[7] );
      rot[1] = rot[1] / mag;
      rot[4] = rot[4] / mag;
      rot[7] = rot[7] / mag;

      mag = sqrt( rot[2]*rot[2] + rot[5]*rot[5] + rot[8]*rot[8] );
      rot[2] = rot[2] / mag;
      rot[5] = rot[5] / mag;
      rot[8] = rot[8] / mag;

   } else */ if( inds[1] == 0 || inds[1] == (grid[0]->sizes[2]-1) || inds[0] == (grid[0]->sizes[0]-1) ) {
      tang1_vec[0]= cos(xi);
      tang1_vec[1] = 0;
      tang1_vec[2] = -1*sin(xi);
      StGermain_VectorCrossProduct(tang2_vec, r_vec, tang1_vec);

      rot[0] = r_vec[0];  rot[1] = tang1_vec[0];  rot[2] = tang2_vec[0];
      rot[3] = r_vec[1];  rot[4] = tang1_vec[1];  rot[5] = tang2_vec[1];
      rot[6] = r_vec[2];  rot[7] = tang1_vec[2];  rot[8] = tang2_vec[2];
   } else {
      ///if (inds[2] == 0 || inds[2] == (grid->sizes[2]-1) )*/ {
      tang1_vec[0] = 0;
      tang1_vec[1] = cos(eta);
      tang1_vec[2] = -sin(eta);
      StGermain_VectorCrossProduct(tang2_vec, r_vec, tang1_vec);

      rot[0] = r_vec[0];  rot[1] = tang2_vec[0];  rot[2] = tang1_vec[0];
      rot[3] = r_vec[1];  rot[4] = tang2_vec[1];  rot[5] = tang1_vec[1];
      rot[6] = r_vec[2];  rot[7] = tang2_vec[2];  rot[8] = tang1_vec[2];
   }


   /*
   double test[9];
   test[0] = rot[0]*rot[0] + rot[1]*rot[1] + rot[2]*rot[2];
   test[1] = rot[0]*rot[3] + rot[1]*rot[4] + rot[2]*rot[5];
   test[2] = rot[0]*rot[6] + rot[1]*rot[7] + rot[2]*rot[8];

   test[3] = rot[3]*rot[0] + rot[4]*rot[1] + rot[5]*rot[2];
   test[4] = rot[3]*rot[3] + rot[4]*rot[4] + rot[5]*rot[5];
   test[5] = rot[3]*rot[6] + rot[4]*rot[7] + rot[5]*rot[8];

   test[6] = rot[6]*rot[0] + rot[7]*rot[1] + rot[8]*rot[2];
   test[7] = rot[6]*rot[3] + rot[7]*rot[4] + rot[8]*rot[5];
   test[8] = rot[6]*rot[6] + rot[7]*rot[7] + rot[8]*rot[8];

   if( fabs(test[0]-1) > 1e-5 ) assert(0);   
   if( fabs(test[4]-1) > 1e-5 ) assert(0);   
   if( fabs(test[8]-1) > 1e-5 ) assert(0);   
   if( test[1] > 1e-5 ) assert(0);
   if( test[2] > 1e-5 ) assert(0);
   if( test[3] > 1e-5 ) assert(0);
   if( test[5] > 1e-5 ) assert(0);
   if( test[6] > 1e-5 ) assert(0);
   if( test[7] > 1e-5 ) assert(0);

   */



}

void Spherical_GetRotationMatrixIJK_FSNodes( Mesh* mesh, int dNodeID, double* rot )
{
   /* Return the rotation matrix of a boundary node which doesn't align with
   * the x-y coordinate system.

definitions:
   x = r * cos(t) * sin(p)
   y = r *          sin(p)
   z = -r * sin(t) * sin(p)



   */

   ProjectionGenerator* gen = NULL;
   Grid **grid = NULL;
   unsigned inds[3];
   double *vert, centre[3], vec[3], dimRes[3], R, x_z_len;

   double r, t, p;
   unsigned gNode;

   centre[0] = centre[1] = centre[2] = 0.0;
   vert = Mesh_GetVertex( mesh, dNodeID );

   // calculate the radius of projection
   r = sqrt( vert[0]*vert[0] + vert[1]*vert[1] + vert[2]*vert[2] );
   assert( r>1e-9 ); // sanity check for zero radius

   x_z_len = sqrt(vert[0]*vert[0] + vert[2]*vert[2]);

   t = atan2( -1*vert[2], vert[0]); // longitude position discretization
   //p = atan2( vert[1], x_z_len );   // latitude position discretization
   p = asin( vert[1]/r );   // latitude position discretization


   rot[0] =      cos(t) * cos(p);
   rot[1] = -r * sin(t) * cos(p);
   rot[2] = -r * cos(t) * sin(p);

   rot[3] =               sin(p);
   rot[4] =  0                  ;
   rot[5] =  r *          cos(p);

   rot[6] = -1 * sin(t) * cos(p);
   rot[7] = -r * cos(t) * cos(p);
   rot[8] =  r * sin(t) * sin(p);
}

#if 0
   phi
   // for 2D case
   if(dim == 2) {
      vec[2] = 0;

      /* assume 2D regular polar mesh, therefore theta(rtp[1]) determines the rotation
      * matrix:
      *   [ cos(theta) -sin(theta) ]
      *   [ sin(theta) cos(theta)  ]

   NOTE: the 2nd column isn't multiplied by the radius of one doesn't have to scale
   between a velocity at the inner and outer radius
  

      */

   } else {

      eta = gen->crdMin[2] + (double)inds[2]*dimRes[2];

      double phi = pow( 1+tan(xi)*tan(xi)+tan(eta)*tan(eta), -0.5 );

      /*
        | dx_dt |     | dx_dxi  dx_deta  dx_dR |  | dxi_dt  |
        | dy_dt |  =  | dy_dxi  dy_deta  dy_dR |  | deta_dt |
        | dz_dt |     | dz_dxi  dz_deta  dz_dR |  | dR_dt   |
         
        here we calculat the coefficient matrix 
      */

      rot[0] = R * phi / (cos(xi)*cos(xi)) * ( 1 - tan(xi)*tan(xi) * phi*phi );
      rot[1] = -R * tan(xi) * tan(eta) * phi*phi*phi / (cos(eta)*cos(eta));
      rot[2] = tan(xi) * phi;

      rot[3] = -R * tan(xi) * tan(eta) * phi*phi*phi / (cos(xi)*cos(xi));
      rot[4] = R * phi / (cos(eta)*cos(eta)) * ( 1 - tan(eta)*tan(eta) * phi*phi );
      rot[5] = tan(eta) * phi;

      rot[6] = -R * phi*phi*phi * tan(xi) / (cos(xi)*cos(xi));
      rot[7] = -R * phi*phi*phi * tan(eta) / (cos(eta)*cos(eta));
      rot[8] = phi;

#if 0
      // just rotate upper surface wall and lower surface
      if( (inds[1] == grid[0]->sizes[1]-1) || inds[1]==0 ) {
      double *vert2=NULL;
      unsigned inds2[3];
      memcpy( inds2, inds, 3*sizeof(unsigned) );
      
         double normal[3];

      if( inds[1] == 0 ) { 
         inds2[1] = 1;
      vert2 = Mesh_GetVertex(mesh, Grid_Project( grid[0], inds2 ));
      StGermain_VectorSubtraction( normal, vert2, vert, 3 );
      } else {
         inds2[1]=inds[1]-1;
      vert2 = Mesh_GetVertex(mesh, Grid_Project( grid[0], inds2 ));
      StGermain_VectorSubtraction( normal, vert, vert2, 3 );
      }

#endif

         double temp[3]={0,0,1};
         double R_1[9], R_2[9];




         memset(R_1, 0, 9*sizeof(double) );
         memset(R_2, 0, 9*sizeof(double) );

#if 0
#endif
         double normal[3];
         //multiply rot by vec to get vector normal to surface
         normal[0] = rot[0]*temp[0]+rot[1]*temp[1]+rot[2]*temp[2];
         normal[1] = rot[3]*temp[0]+rot[4]*temp[1]+rot[5]*temp[2];
         normal[2] = rot[6]*temp[0]+rot[7]*temp[1]+rot[8]*temp[2];

         double angle, angle2;
         double mag = sqrt(normal[0]*normal[0]+normal[1]*normal[1]+normal[2]*normal[2]);
         double mag_x_y=sqrt(normal[0]*normal[0]+normal[1]*normal[1]);
         if( mag_x_y < 1e-15 ) {
            angle=0;//mag_x_y=mag;
         }
         else {
            angle=acos(normal[1]/mag_x_y);
         }
         if( normal[0] > 0 ) { angle = -1*angle; }
         R_1[0]=cos(angle); R_1[1]=-sin(angle); R_1[2]=0;
         R_1[3]=sin(angle); R_1[4]=cos(angle); R_1[5]=0;
         R_1[6]=0; R_1[7]=0; R_1[8]=1;

         angle2=acos(mag_x_y/mag);
         R_2[0]=1;
         R_2[4]=cos(angle2); R_2[5]=-sin(angle2);
         R_2[7]=sin(angle2); R_2[8]=cos(angle2);

         blasMatrixMult( R_1, R_2, 3, 3, 3, rot );

         // final rotation of local x along east-west axis
         memcpy(R_2, rot, 9*sizeof(double) );
         double x_unit[3] = {1,0,0}; //true east
         double rot_x_vector[3]; 

         // build the current x basis vector for our rotate coord system
         rot_x_vector[0] = rot[0]; 
         rot_x_vector[1] = rot[3]; 
         rot_x_vector[2] = rot[6]; 

         // find difference btw true x and current x basis
         angle2 = StGermain_AngleBetweenVectors( rot_x_vector, x_unit, 3 );

         // reverse angles if xi > 0, condition due to previous rotations
         if( xi>0 )
            angle2 = -1*angle2;

         memset(R_1, 0, 9*sizeof(double) );
         R_1[0]=cos(angle2); R_1[1]=0; R_1[2]=-sin(angle2);
         R_1[3]=0; R_1[4]=1; R_1[5]=0;
         R_1[6]=sin(angle2); R_1[7]=0; R_1[8]=cos(angle2);

         blasMatrixMult( R_2, R_1, 3, 3, 3, rot );

      }
#if 0
      if( inds[1] == 0 ) {
         //east wall
      } else if ( inds[0] == gen->sizes ) {

      }
#endif
#if 0
      double w;
      // intermediate var
      w = 1 + tan(xi)*tan(xi) + tan(eta)*tan(eta);

      rot[0] =  R/(cos(xi)*cos(xi))*pow(w,-0.5) * (1-tan(xi)*tan(xi)/w);
      rot[1] =  -R*tan(xi)*tan(eta)/(cos(eta)*cos(eta)) * pow(w, -1.5);
      rot[2] = tan(xi) * pow(w,-0.5);

      rot[3] =  -R*tan(xi)*tan(eta)/(cos(xi)*cos(xi)) * pow(w,-1.5);
      rot[4] =  R/(cos(eta)*cos(eta)) * pow(w,-0.5) * (1-tan(eta)*tan(eta)/w);
      rot[5] =  tan(eta) * pow(w,-0.5);

      rot[6] = -R*tan(xi)/(cos(xi)*cos(xi)) * pow(w,-1.5);
      rot[7] = -R*tan(eta)/(cos(eta)*cos(eta)) * pow(w,-1.5);
      rot[8] = pow(w,-0.5);
#endif

#endif

void Spherical_GetRotationMatrixIJK_SphericalNodes( Mesh* mesh, int dNodeID, double* rot )
{
   /* 
      Purpose: Return the rotation matrix of a boundary node to transform its velocity from
               the x-y coordinate system to radius-theta-phi
    */

   double *vert, centre[3], vec[3];
   unsigned dim = mesh->generator->nDims;

   assert( rot );

   // TODO: simple centre definition for now, should make more flexible
   centre[0] = 0.0;
   centre[1] = 0.0;
   centre[2] = 0.0;

   vert = Mesh_GetVertex( mesh, dNodeID );

   // calculate the distance between n1 and n2
   StGermain_VectorSubtraction( vec, vert, centre, dim );

   // for 2D case
   if(dim == 2) {
      double r, t;
      unsigned int inds[3];
      unsigned gNode;

      // get the vert grid and generator for quick lookup
      Grid **grid = (Grid** )Mesh_GetExtension(mesh,Grid**,mesh->vertGridId);
      SphericalGenerator *sg=(SphericalGenerator*)mesh->generator;

      gNode = Mesh_DomainToGlobal( mesh, MT_VERTEX, dNodeID );
      Grid_Lift( grid[0], gNode, inds );

      r = sg->crdMin[0] + sg->sph_res[0] * (double)inds[0];
      t = (M_PI/180.0) * (sg->crdMin[1] + sg->sph_res[1]*(double)inds[1]); // longitude position discretization

      /* assume 2D regular polar mesh, therefore theta(rtp[1]) determines the rotation
      * matrix:
      *   [ cos(theta) -sin(theta) ]
      *   [ sin(theta) cos(theta)  ]
      */
      rot[0] =      cos(t);
      rot[1] = -r * sin(t);
      rot[2] =      sin(t);
      rot[3] =  r * cos(t);
   } else {
      double r, t, p;
      unsigned int inds[3];
      unsigned gNode;

      // get the vert grid and generator for quick lookup
      Grid **grid = (Grid** )Mesh_GetExtension(mesh,Grid**,mesh->vertGridId);
      SphericalGenerator *sg=(SphericalGenerator*)mesh->generator;

      gNode = Mesh_DomainToGlobal( mesh, MT_VERTEX, dNodeID );
      Grid_Lift( grid[0], gNode, inds );

      r = sg->crdMin[0] + sg->sph_res[0] * (double)inds[0];
      t = (M_PI/180.0) * (sg->crdMin[1] + sg->sph_res[1]*(double)inds[1]); // longitude position discretization
      p = (M_PI/180.0) * (sg->crdMin[2] + sg->sph_res[2] *(double)inds[2]);   // latitude position discretization


      rot[0] =      cos(t) * cos(p);
      rot[1] = -r * sin(t) * cos(p);
      rot[2] = -r * cos(t) * sin(p);

      rot[3] =               sin(p);
      rot[4] =  0                  ;
      rot[5] =  r *          cos(p);

      rot[6] = -1 * sin(t) * cos(p);
      rot[7] = -r * cos(t) * cos(p);
      rot[8] =  r * sin(t) * sin(p);

#if 0
      double U[9], vT[9], sValues[3];
      double *work, optwork;
      int order, lda, ldu, ldvt, LWORK, IWORK, INFO;
      char a='A';
      char trans='T', n='N';
      double C[9], D[9], one, zero;
      one=1; zero=0;
      order=ldu=lda=ldvt=3;

      Spherical_XYZ2RTP3D( vec, rtp );
      r = rtp[0]; // radius
      t = rtp[1]; // angle theta
      p = rtp[2]; // angle phi


      rot[0] =      sin(t) * cos(p);
      rot[1] =  r * cos(t) * cos(p);
      rot[2] = -r * sin(t) * sin(p);

      rot[3] =               sin(p);
      rot[4] =  0;
      rot[5] =  r *          cos(p);

      rot[6] =      cos(t) * cos(p);
      rot[7] = -r * sin(t) * cos(p);
      rot[8] = -r * cos(t) * sin(p);


      // first find optimal work space
      LWORK=-1;
      dgesdd_( &a, &order, &order, rot, &lda, 
            sValues, U, &ldu, vT, &ldvt, 
            &optwork, &LWORK, IWORK, &INFO );

      // allocate work space
      LWORK=(int)optwork;
      work = (double*)malloc( LWORK*sizeof(double) );

      // SVD compute
      dgesdd_( &a, &order, &order, rot, &lda, 
            sValues, U, &ldu, vT, &ldvt, 
            work, &LWORK, IWORK, &INFO );

      if( INFO != 0 ) {
        assert(0);
      }

      dgemm_(&n,&trans,&order, &order, &order, &one, vT, &order, U, &order, &zero, rot, &order );
      dgemm_(&trans,&n,&order, &order, &order, &one, rot, &order, rot, &order, &zero, D, &order );

      // check if D is the identity matrix
      // algorithm is look at every element. If 
      int i,j;
      for(i=0;i<order;i++) {
         for(j=0;j<order;j++) {
            // if element is non-zero
            if( fabs(D[i*order+j]) > 1e-5 ) {
               // if on diagonal and one, that's cool
               if( i==j && fabs(D[i*order+j]-1)<1e-5 ) { continue; }
               else {
                  // in not... bad 
                  printf("strange\n");
                  assert(0);
               }
            }
         }
      }

      free(work);
#endif

   }

}


void Spherical_FeVariable_NonAABCsCalibration( void *_self )
{
   /*@

   	This function rotates the boundary node values of an feVariable into xyz coordinates

   	Assume - the rotation is defined by Spherical_GetRotationMatrixIJK()

   	@*/

   FeVariable *fe = (FeVariable*)_self;
   FeMesh *mesh = fe->feMesh;
   Node_LocalIndex     node_I;
   Node_LocalIndex     nodeCount;
   double              vec[3], vec_xy[3], rot[9];

   assert( mesh );

   nodeCount = FeMesh_GetNodeDomainSize( fe->feMesh );

   for ( node_I = 0 ; node_I < nodeCount ; node_I++ )
   {
      /* if not a velocity boundary nodes skip */
      if( !IndexSet_IsMember( mesh->bndNodeSet, node_I ) )
      {
         continue;
      }

      /* get the veloccity defined in tn */
      FeVariable_GetValueAtNode( fe, node_I, vec );
      /* get rotation matrix for the boundary node */
      Spherical_Get_RotationMatrixIJK( (Mesh*)fe->feMesh, node_I, rot );

      /** vec_xy = R * vec_tn */
      if( mesh->generator->nDims == 2 ) {
         vec_xy[0] = rot[0]*vec[0] + rot[1]*vec[1];
         vec_xy[1] = rot[2]*vec[0] + rot[3]*vec[1];
      } else {
         vec_xy[0] = rot[0]*vec[0] + rot[1]*vec[1] + rot[2]*vec[2];
         vec_xy[1] = rot[3]*vec[0] + rot[4]*vec[1] + rot[5]*vec[2];
         vec_xy[2] = rot[6]*vec[0] + rot[7]*vec[1] + rot[8]*vec[2];
      }
      /* set velocity */
      FeVariable_SetValueAtNode( fe, node_I, vec_xy );
   }
}

void Spherical_VectorRTP2XYZ( double *Q, double *xyz, int dim, double* v ) {
   /*@
   	Purpose - Convert spherical vector, Q, at position xyz, into cartesian vector, v, assuming
		  center of spherical coordinates is at x=0, y=0

		     | dQ_dx |     = | dQ_dr  dQ_dtheta | | dr_dx |
		     | dQ_dy |       | dQ_dr  dQ_dtheta | | dtheta_dy |
   @*/

   double radius;
   double norm_cart, norm_rtp;
   
   if( dim == 2 ) {
      radius = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] );

      v[0] = ( Q[0]*xyz[0] - Q[1]*xyz[1] ) / radius;
      v[1] = ( Q[0]*xyz[1] + Q[1]*xyz[0] ) / radius;

      norm_cart = sqrt( Q[0]*Q[0] + Q[1]*Q[1] );
      norm_rtp  = sqrt( v[0]*v[0] + v[1]*v[1] );


   } else {assert(0);}

   if( fabs(norm_cart - norm_rtp) > 1e-5 ) {
      printf("\nFucked up in %s\n\n", __func__ );
      assert(0);
   }
}

void Spherical_VectorXYZ2RTP( double *_v, double *xyz, int dim, double* v )
{
   /*@
   	Purpose - Convert cartesian vector, _v, at position xyz, into spherical vector, v, assuming
		  center of spherical coordinates is at x=0, y=0

		     | v_r |     = | dr_dx  dr_dy | | dx_dt |
		     | v_T |       | dT_dx  dT_dy | | dy_dt |

		  here v_r & v_T are the spherical velocity components radius and theta respectively
		  dx_dt & dy_dt are the cartesian velocity components
   @*/

   double radius;
   double norm_cart, norm_rtp;
   
   if( dim == 2 ) {
      radius = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] );

      v[0] = ( _v[0]*xyz[0] + _v[1]*xyz[1] ) / radius;
      v[1] = (-_v[0]*xyz[1] + _v[1]*xyz[0] ) / radius;


      norm_cart = sqrt( _v[0]*_v[0]+_v[1]*_v[1] );
      norm_rtp = sqrt( v[0]*v[0] + v[1]*v[1] );
   } 
   else {

      double b1, b2, radius;
      radius = sqrt( xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
      b2 = xyz[0]*xyz[0] + xyz[2]*xyz[2];

      if( fabs(b2)<1e-13 ) assert(0);
      b1 = xyz[1] / sqrt(b2);

      v[0] = (_v[0]*xyz[0] + _v[1]*xyz[1] + _v[2]*xyz[2] ) / radius;

      v[1] = (_v[0] * -xyz[2] + xyz[0]*_v[2]) / (b2);

      v[2] = ( -xyz[0]*xyz[1]*pow(b2, -1.5) * _v[0] +
             pow(b2, -0.5) * _v[1] +
             -xyz[1]*xyz[2] * pow(b2,-1.5) ) / (1 + b1*b1 );

      norm_cart = sqrt( _v[0]*_v[0]+_v[1]*_v[1]+_v[2]*_v[2] );
      norm_rtp = sqrt( v[0]*v[0] + v[1]*v[1]+v[2]*v[2] );
   }

   if( fabs(norm_cart - norm_rtp) > 1e-5 ) {
      printf("\nFucked up in %s\n\n", __func__ );
      assert(0);
   }



}


void Spherical_XYZ2RTP2D( double *xyz, double* rtp )
{
   /*@
   	Purpose - Convert cartesian coordinates in spherical coordinates assuming
   	center of spherical coordinates is at x=0, y=0

        Return rtp vector is [radius, theta], where theta is between [-pi, pi]
   @*/

   // sort out radius
   rtp[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] );

   rtp[1] = atan2( xyz[1], xyz[0] );
}

void Spherical_XYZ2RTP3D( double *xyz, double* rtp )
{
   /*@
     Don't know if implementation is good
   	Purpose - Convert cartesian coordinates in spherical coordinates assuming
   	center of spherical coordinates is at x=0, y=0, z=0

      assumes:
       x = r * cos(t) * cos(p)
       y = r *          sin(p)
       z =-r * sin(t) * cos(p)
      where r is the radius, t is the longitude and p is the latitude

   @*/
   // sort out radius
   rtp[0] = sqrt(xyz[0]*xyz[0] + xyz[1]*xyz[1] + xyz[2]*xyz[2] );
   rtp[1] = atan2( -1*xyz[2], xyz[0]);
   rtp[2] = asin( xyz[1]/rtp[0] );

#if 0
   if( fabs(xyz[0]) < DBL_MIN ) {
      if( xyz[2] > 0 ) rtp[1] = -M_PI/2;
      else rtp[1] = M_PI/2;
   } else
      rtp[1] = atan( -1*xyz[2]/xyz[0] );

   /* the new projection mesh stuff. Was only concerntrating on the sixth of a mesh */
   // sort out phi
   rtp[2] = asin(xyz[2]/rtp[0]);

   assert( rtp[2] <= M_PI/2.0 );
   assert( rtp[2] >= -M_PI/2.0 );

   // sort out theta
   rtp[1] = asin( xyz[1]/(rtp[0]*cos(rtp[2]) ) );

  /* Feb 28th 2013
      old implementation before the new projection sphere mesh stuff
      */

   // sort out theta
   if (fabs(xyz[0]) < DBL_MIN ) {
      if( xyz[1] > DBL_MIN )	       rtp[1] = M_PI/2;
      else if ( xyz[1] < -1*DBL_MIN )  rtp[1] = 3*M_PI/2;
      else			       rtp[1] = 0;
   }
   else
   {
      angle = atan( xyz[1]/xyz[0] );
      if( xyz[0] > 0 ) {
	 if( xyz[1] > 0 )  rtp[1] = angle;	      // north-east quad
	 else		   rtp[1] = 2*M_PI - angle;   // south-east quad
      }
      else
      {
	 if( xyz[1] > 0 )  rtp[1] = M_PI - angle;     // north-west 
	 else		   rtp[1] = M_PI + angle;     // south-west
      }
   }

   // sort out phi

   if (fabs(xyz[0]) < DBL_MIN ) {
      if( xyz[2] > DBL_MIN )		  rtp[2] = M_PI/2;
      else if ( xyz[2] < -1*DBL_MIN )	  rtp[2] = 3*M_PI/2;
      else				  rtp[2] = 0;
   }
   else
   {
      angle = atan( xyz[2]/xyz[0] );
      if( xyz[0] > 0 ) {
	 if( xyz[2] > 0 )  rtp[2] = angle;	      // north-east quad
	 else		   rtp[2] = 2*M_PI - angle;   // south-east quad
      }
      else
      {
	 if( xyz[2] > 0 )  rtp[2] = M_PI - angle;     // north-west 
	 else		   rtp[2] = M_PI + angle;     // south-west
      }
   }
#endif
   
}

void Spherical_RTP2XYZ( double *rtp, double* xyz )
{
   /*@
   	Purpose - Convert spherical polar coordinates in cartesiancoordinates
   @*/

   xyz[0] = rtp[0] * cos( rtp[1] );
   xyz[1] = rtp[0] * sin( rtp[1] );
}

void Cylindrical_WithPerturbation(Node_LocalIndex node_lI,Variable_Index var_I,void *_context,void* _data, void* _result) {
   /* This is the IC condition used for the Cylindrical Benchmark.
      It prescribes a temperature IC that is linear in depth with a periodic perturbation 
   */
   
   UnderworldContext* context = (UnderworldContext*)_context;
   FeVariable* feVariable = (FeVariable*)FieldVariable_Register_GetByName( context->fieldVariable_Register, "TemperatureField" );
   FeMesh *feMesh = feVariable->feMesh;
   double rtp[3];
   double* result = (double*) _result;

   unsigned nDims = Mesh_GetDimSize( feMesh );

   double *coord = Mesh_GetVertex( feMesh, node_lI );

   Spherical_XYZ2RTP2D( coord, rtp );

   *result=rtp[1];


#if 0
   Mesh_GetGlobalCoordRange( feMesh, min, max );

   topLayerCoord = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerCoord", max[vertaxis] );
   bottomLayerCoord = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerCoord", min[vertaxis] );

   topLayerBC = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_TopLayerBC", 0.0 );
   bottomLayerBC = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_BottomLayerBC", 1.0 );
   scaleFactor = bottomLayerBC - topLayerBC;
   perturbationAmplitude = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_PerturbationAmplitude", 0.1 );
   /* Note: these are both multiplied by pi, so wavenumber = 1 means the perturbation goes from 0 to pi, which is
    * half a full sin or cos cycle. Wavenumber = 3 means the range is 0 -> 3pi, or 1 and a half full cycles. */
   horizontalWaveNumber = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_HorizontalWaveNumber", 1.0 );
   verticalWaveNumber = Dictionary_GetDouble_WithScopeDefault( dictionary, (Dictionary_Entry_Key)"SinusoidalTempIC_VerticalWaveNumber", 1.0 );


   /* if node is outside IC shape set to 0 temperature */
   if( coord[vertaxis] > topLayerCoord || coord[vertaxis] < bottomLayerCoord ) {
      *result = 0; return ;
   }

   /* make coord relative to box bottom left corner, then scale from 0 to 1 between box min & max */
   relScaledCoord[I_AXIS]   = (coord[0] - min[0]) / (max[0] - min[0]);
   relScaledCoord[vertaxis] = (coord[vertaxis] - bottomLayerCoord) / (topLayerCoord - bottomLayerCoord);


   /* Note: ok to use the 1.0 below since we've already scaled the coord to somewhere between 0 to 1 */
   pertCoeff = cos( horizontalWaveNumber * M_PI * relScaledCoord[ I_AXIS ] ) * sin( verticalWaveNumber * M_PI * relScaledCoord[ vertaxis ] ) ;
   *result = topLayerBC + scaleFactor * (1.0 - relScaledCoord[ vertaxis ]) + perturbationAmplitude * pertCoeff ;
#endif

}
