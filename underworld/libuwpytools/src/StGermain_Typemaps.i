/* -*- C -*-  (not really, but good for syntax highlighting) */

%{
/* Includes the header in the wrapper code */
#include <StGermain/StGermain.h>

%}

/* we include this rename for Variable.h which has a member 'as' which conflicts with the 'as' python builtin */
%rename(atType) as;

%typemap(out) Stg_Class* {

   if ($1==0)
      return NULL; 

   PyObject *obj=NULL;

   Stg_Class* component = (Stg_Class*) $1;
   char* typename = (char*)malloc(sizeof(char) * (strlen(component->type)+3));
   sprintf(typename, "%s *", component->type);
   swig_type_info *swigTypeInfo = SWIG_Python_TypeQuery(typename);
   if(swigTypeInfo)
      obj = SWIG_NewPointerObj((void*)$1, swigTypeInfo, 0);
   else
      obj = SWIG_NewPointerObj((void*)$1, $1_descriptor, 0);
   free(typename);

   $result = obj;
}

%typemap(out) Stg_Object* {

   if ($1==0)
      return NULL; 

   PyObject *obj=NULL;

   Stg_Class* component = (Stg_Class*) $1;
   char* typename = (char*)malloc(sizeof(char) * (strlen(component->type)+3));
   sprintf(typename, "%s *", component->type);
   swig_type_info *swigTypeInfo = SWIG_Python_TypeQuery(typename);
   if(swigTypeInfo)
      obj = SWIG_NewPointerObj((void*)$1, swigTypeInfo, 0);
   else
      obj = SWIG_NewPointerObj((void*)$1, $1_descriptor, 0);
   free(typename);

   $result = obj;
}



%typemap(out) Stg_Component* {

   if ($1==0)
      return NULL; 

   PyObject *obj=NULL;

   Stg_Class* component = (Stg_Class*) $1;
   char* typename = (char*)malloc(sizeof(char) * (strlen(component->type)+3));
   sprintf(typename, "%s *", component->type);
   swig_type_info *swigTypeInfo = SWIG_Python_TypeQuery(typename);
   if(swigTypeInfo)
      obj = SWIG_NewPointerObj((void*)$1, swigTypeInfo, 0);
   else
      obj = SWIG_NewPointerObj((void*)$1, $1_descriptor, 0);
   free(typename);

   $result = obj;
}


