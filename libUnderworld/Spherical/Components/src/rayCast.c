/* all code, expect for interface function RayCast_IsPointInElement,
   has been taken from

   http://rosettacode.org/wiki/Ray-casting_algorithm
   under the GNU Free Documentation License 1.2
*/


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <mpi.h>
#include <float.h>
#include <math.h>

#include <StGermain/StGermain.h>
#include <StgDomain/StgDomain.h>
#include <StgFEM/StgFEM.h>
#include <PICellerator/PICellerator.h>

#include "Components.h"


 
#define BIN_V(op, xx, yy) vec v##op(vec a,vec b){vec c;c.x=xx;c.y=yy;return c;}
#define BIN_S(op, r) double v##op(vec a, vec b){ return r; }
BIN_V(sub, a.x - b.x, a.y - b.y);
BIN_V(add, a.x + b.x, a.y + b.y);
BIN_S(dot, a.x * b.x + a.y * b.y);
BIN_S(cross, a.x * b.y - a.y * b.x);
 
/* return a + s * b */
vec vmadd(vec a, double s, vec b)
{
	vec c;
	c.x = a.x + s * b.x;
	c.y = a.y + s * b.y;
	return c;
}
 
/* check if x0->x1 edge crosses y0->y1 edge. dx = x1 - x0, dy = y1 - y0, then
   solve  x0 + a * dx == y0 + b * dy with a, b in real
   cross both sides with dx, then: (remember, cross product is a scalar)
	x0 X dx = y0 X dx + b * (dy X dx)
   similarly,
	x0 X dy + a * (dx X dy) == y0 X dy
   there is an intersection iff 0 <= a <= 1 and 0 <= b <= 1
 
   returns: 1 for intersect, -1 for not, 0 for hard to say (if the intersect
   point is too close to y0 or y1)
*/
int intersect(vec x0, vec x1, vec y0, vec y1, double tol, vec *sect)
{
	vec dx = vsub(x1, x0), dy = vsub(y1, y0);
	double d = vcross(dy, dx), a;
	if (!d) return 0; /* edges are parallel */
 
	a = (vcross(x0, dx) - vcross(y0, dx)) / d;
	if (sect)
		*sect = vmadd(y0, a, dy);
 
	if (a < -tol || a > 1 + tol) return -1;
	if (a < tol || a > 1 - tol) return 0;
 
	a = (vcross(x0, dy) - vcross(y0, dy)) / d;
	if (a < 0 || a > 1) return -1;
 
	return 1;
}
 
/* distance between x and nearest point on y0->y1 segment.  if the point
   lies outside the segment, returns infinity */
double dist(vec x, vec y0, vec y1, double tol)
{
	vec dy = vsub(y1, y0);
	vec x1, s;
	int r;
 
	x1.x = x.x + dy.y; x1.y = x.y - dy.x;
	r = intersect(x, x1, y0, y1, tol, &s);
	if (r == -1) return HUGE_VAL;
	s = vsub(s, x);
	return sqrt(vdot(s, s));
}
 
#define for_v(i, z, p) for(i = 0, z = p->v; i < p->n; i++, z++)
/* returns 1 for inside, -1 for outside, 0 for on edge */
int inside(vec v, polygon p, double tol)
{
	/* should assert p->n > 1 */
	int i, k, crosses;
	vec *pv;
	double min_x, max_x, min_y, max_y;
 
	for (i = 0; i < p->n; i++) {
		k = i + 1 % p->n;
		min_x = dist(v, p->v[i], p->v[k], tol);
		if (min_x < tol) return 0;
	}
 
	min_x = max_x = p->v[0].x;
	min_y = max_y = p->v[1].y;
 
	/* calculate extent of polygon */
	for_v(i, pv, p) {
		if (pv->x > max_x) max_x = pv->x;
		if (pv->x < min_x) min_x = pv->x;
		if (pv->y > max_y) max_y = pv->y;
		if (pv->y < min_y) min_y = pv->y;
	}
	if (v.x < min_x || v.x > max_x || v.y < min_y || v.y > max_y)
		return -1;
 
	max_x -= min_x; max_x *= 2;
	max_y -= min_y; max_y *= 2;
	max_x += max_y;
 
	vec e;
	while (1) {
		crosses = 0;
		/* pick a rand point far enough to be outside polygon */
		e.x = v.x + (1 + rand() / (RAND_MAX + 1.)) * max_x;
		e.y = v.y + (1 + rand() / (RAND_MAX + 1.)) * max_x;
 
		for (i = 0; i < p->n; i++) {
			k = (i + 1) % p->n;
			k = intersect(v, e, p->v[i], p->v[k], tol, 0);
 
			/* picked a bad point, ray got too close to vertex.
			   re-pick */
			if (!k) break;
 
			if (k == 1) crosses++;
		}
		if (i == p->n) break;
	}
	return (crosses & 1) ? 1 : -1;
}

Bool RayCast_IsPointInFace( FeMesh* mesh, unsigned elInd, const double* coord ){
   IArray *inc=mesh->inc; 
   double *vcoord; // pointers to vertex coords
   double lengthScale;
   vec vecList[4];
   int *nodes, nNodes, nDims;
   int i;

   // algorithm inputs
   polygon_t poly;
   vec myVec;

   // get nodes of element - now call them vertexes
   nDims = Mesh_GetDimSize( mesh );
   Mesh_GetIncidence( mesh, MT_FACE, elInd, MT_VERTEX, inc );
   nNodes = IArray_GetSize( inc );
   nodes = IArray_GetPtr( inc );

   if( fabs(nodes[0]-nodes[1]) > 1 ) {
      printf("Say what! node list is weird\n");
   }
   /* build vertex list, not we ASSUME a 2D element has
      an connectivity of node 0->1->3->2
   */

   i = 0;
   vcoord = Mesh_GetVertex( mesh, nodes[0] );
   vecList[i].x = vcoord[0];
   vecList[i].y = vcoord[1];

   i = 1;
   vcoord = Mesh_GetVertex( mesh, nodes[1] );
   vecList[i].x = vcoord[0];
   vecList[i].y = vcoord[1];

   i = 2;
   vcoord = Mesh_GetVertex( mesh, nodes[3] ); // HERE IS THE IMPORTANT SWITCH
   vecList[i].x = vcoord[0];
   vecList[i].y = vcoord[1];

   i = 3;
   vcoord = Mesh_GetVertex( mesh, nodes[2] ); // HERE IS THE IMPORTANT SWITCH
   vecList[i].x = vcoord[0];
   vecList[i].y = vcoord[1];

   // measure length scale
   lengthScale = fabs(vecList[0].x-vecList[2].x);

   poly.n = 4;
   poly.v = vecList;

   myVec.x = coord[0];
   myVec.y = coord[1];

   i = inside( myVec, &poly, 1e-8*lengthScale );
   if( i < 0 )
      return False;
   else 
      return True;
}

Bool RayCast_IsPointIn3DElement( FeMesh* mesh, unsigned elInd, const double* coord ){
   /* 
     If we have a 3D element we perform 6 raycast checks
     one on each face of the element (irregular hexahedron).
     If the point is within each face, then it must be inside the element
   */

   IArray *inc=mesh->inc; 
   double *vcoord; // pointers to vertex coords
   double lengthScale;
   vec vecList[4];
   int *nodes, nNodes, nDims;
   int i,nFaces,*faces, f_i;

   IArray* face_inc=IArray_New();
   // algorithm inputs
   polygon_t poly;
   vec myVec;

   // get nodes of element - now call them vertexes
   nDims = Mesh_GetDimSize( mesh );
   Mesh_GetIncidence( mesh, Mesh_GetDimSize(mesh), elInd, MT_FACE, face_inc );
   nFaces = IArray_GetSize( face_inc );
   faces = IArray_GetPtr( face_inc );

   for( f_i=0;f_i<nFaces;f_i++ ) {
      // get the nodes in face
      Mesh_GetIncidence(mesh, MT_FACE, faces[f_i], MT_VERTEX, inc ); 
      
   }

   if( fabs(nodes[0]-nodes[1]) > 1 ) {
      printf("Say what! node list is weird\n");
   }
   /* build vertex list, not we ASSUME a 2D element has
      an connectivity of node 0->1->3->2
   */
   Stg_Class_Delete( face_inc );

}
