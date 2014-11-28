typedef struct { double x, y; } vec;
typedef struct { int n; vec* v; } polygon_t, *polygon;

/* returns 1 for inside, -1 for outside, 0 for on edge */
int inside(vec v, polygon p, double tol);

/* 
FUNC: Underworld interface for a 2D ray casting algorithm to check 
   if a global coordinate is inside an element

RETURNS: 
   True ... coord is inside or on the edge of element elInd
   False ... coord is outside element elInd

INPUT:
   *mesh ... a mesh pointer
   elInd ... the global index of the element to test
   *coord ... coordinates to test

ASSUMES: 
   element node ordering in 2D is 0->1->3->2
*/
Bool RayCast_IsPointInElement(FeMesh* mesh, unsigned elInd, const double* coord );
