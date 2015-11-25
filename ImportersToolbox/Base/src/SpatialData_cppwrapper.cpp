/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c**
** Copyright (C), 2003, Victorian Partnership for Advanced Computing (VPAC) Ltd, 110 Victoria Street, Melbourne, 3053, Australia.
**
** Authors:
** Stevan M. Quenette, Senior Software Engineer, VPAC. (steve@vpac.org)
** Patrick D. Sunter, Software Engineer, VPAC. (pds@vpac.org)
** Luke J. Hodkinson, Computational Engineer, VPAC. (lhodkins@vpac.org)
** Siew-Ching Tan, Software Engineer, VPAC. (siew@vpac.org)
** Alan H. Lo, Computational Engineer, VPAC. (alan@vpac.org)
** Raquibul Hassan, Computational Engineer, VPAC. (raq@vpac.org)
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Lesser General Public
**  License as published by the Free Software Foundation; either
**  version 2.1 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Lesser General Public License for more details.
**
**  You should have received a copy of the GNU Lesser General Public
**  License along with this library; if not, write to the Free Software
**  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
**
**~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
#ifdef HAVE_SPATIALDATA

#include "SpatialData_cppwrapper.hh"

#include "spatialdata/spatialdb/SimpleDB.hh"
#include "spatialdata/spatialdb/SimpleIO.hh"
#include "spatialdata/spatialdb/SimpleIOAscii.hh"

#include "spatialdata/geocoords/CoordSys.hh" 
#include "spatialdata/geocoords/CSCart.hh"

#include "spatialdata/spatialdb/SimpleDBData.hh"
#include <stdio.h>
#include <string.h>

typedef struct  {
   spatialdata::spatialdb::SimpleDB* simpleDB;
   spatialdata::spatialdb::SimpleIOAscii* dbIO;
   spatialdata::geocoords::CSCart* csCart;
   int numVals;
} SDdata;

extern "C" int _SpatialData_cppwrapper_Build( void** data, char* filename, char* approxType, char* fieldName ) {
   /** init the db */
   SDdata* dataSD = new SDdata;
   dataSD->simpleDB = new spatialdata::spatialdb::SimpleDB;
   dataSD->dbIO     = new spatialdata::spatialdb::SimpleIOAscii;
   dataSD->dbIO->filename(filename);
   dataSD->simpleDB->ioHandler(dataSD->dbIO);

   dataSD->simpleDB->open();

   if(!strcasecmp("linear", approxType))
      dataSD->simpleDB->queryType(spatialdata::spatialdb::SimpleDB::LINEAR);
   else
      dataSD->simpleDB->queryType(spatialdata::spatialdb::SimpleDB::NEAREST);

   const char* names[] = {fieldName};
   dataSD->numVals = 1;
 
   dataSD->simpleDB->queryVals(names, dataSD->numVals);

   dataSD->csCart = new spatialdata::geocoords::CSCart;
   dataSD->csCart->initialize();
   *data = (void*)dataSD;

   return 1;
}

extern "C" int _SpatialData_cppwrapper_Destroy( void** data ) {
   SDdata* dataSD = (SDdata*) *data;
   // close db 
   dataSD->simpleDB->close();
   delete dataSD->dbIO;
   delete dataSD->simpleDB;
   delete dataSD;

   return 1;
}


extern "C" int _SpatialData_cppwrapper_evaluate( void** data, double* queryLoc, double* value ) {
   SDdata* dataSD = (SDdata*) *data;
   //const int errFlags[] = { 0 };
   const int spaceDim = 3;
   
   return dataSD->simpleDB->query(value, dataSD->numVals, queryLoc, spaceDim,
                dataSD->csCart);

}




#endif


