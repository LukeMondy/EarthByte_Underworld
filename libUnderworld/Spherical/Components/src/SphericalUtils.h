/* See SphericalUtils.c file for function documentation */

//void Spherical_GetRotationMatrixIJK( Mesh* mesh, int dNodeID, double* rot );

void Spherical_GetRotationMatrixIJK_RSNodes( Mesh* mesh, int dNodeID, double* rot );
void Spherical_GetRotationMatrixIJK_FSNodes( Mesh* mesh, int dNodeID, double* rot );
void Spherical_GetRotationMatrixIJK_SphericalNodes( Mesh* mesh, int dNodeID, double* rot );
void Spherical_FeVariable_NonAABCsCalibration( void *fe );

void Spherical_RTP2XYZ( double *rtp, double* xyz );

void Spherical_XYZ2RTP( double *xyz, double* rtp );

/* Returns spherical vector v. Converts cartesian vector, _v, at position xyz assuming center of spherical coordinates is at x=0, y=0 */
void Spherical_VectorXYZ2RTP( double *_v, double *xyz, int dim, double* v );

void Spherical_XYZ2RTP2D( double *xyz, double* rtp );
void Spherical_XYZ2RTP3D( double *xyz, double* rtp );
void (*Spherical_Get_RotationMatrixIJK)(Mesh* mesh, int dNodeID, double* rot);

void Cylindrical_WithPerturbation(Node_LocalIndex node_lI,Variable_Index var_I,void *_context,void* _data, void* _result);
