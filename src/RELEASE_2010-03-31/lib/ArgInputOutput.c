#include "GRACEiolib.h"
#include "GRACEio_utils.h"
#include <string.h>

static int8_t SccsId[] = "$Id: ArgInputOutput.c,v 1.8 2009/11/16 21:21:30 glk Exp $";

/* ************************************************************************** */
int32_t FindOpt (int32_t argc, int8_t *argv[], int8_t *str)

{
  int32_t        i, res;

  res = -1;
  for (i=1; i< argc; i++)
     if (strcmp(argv[i], str) == 0)
     res= i;

  return(res);
}

/* ************************************************************************** */
int32_t GetLongArgv(int32_t argc,int8_t *argv[], int8_t *str)

{

 int32_t        idum;

 int32_t       LongVal;

 idum = FindOpt(argc,argv, str);

 if (idum != -1 && idum <= argc-2)
   {
    sscanf(argv[idum + 1], "%d", &LongVal);
    return LongVal;
   }
 else
   {
    fprintf(stderr,"\n %s must be specified ! \n\n",str);
    exit(1);
   }
}

void Get2LongArgv(int32_t argc,int8_t *argv[], int8_t *str, int32_t *quan1, int32_t *quan2)

{

 int32_t        idum;

 idum = FindOpt(argc,argv, str);

 if (idum != -1 && idum <= argc-3)
   {
    sscanf(argv[idum + 1], "%d", quan1);
    sscanf(argv[idum + 2], "%d", quan2);
    return ;
   }
 else
   {
    fprintf(stderr,"\n %s must be specified ! \n\n",str);
    exit(1);
   }
}
/* ************************************************************************** */
int32_t GetIntArgv(int32_t argc,int8_t *argv[], int8_t *str)

{

 int32_t        idum;

 int32_t        IntVal;

 idum = FindOpt(argc,argv, str);

 if (idum != -1 && idum <= argc-2)
   {
    sscanf(argv[idum + 1], "%d", &IntVal);
    return IntVal;
   }
 else
   {
    fprintf(stderr,"\n %s must be specified ! \n\n",str);
    exit(1);
   }
}
/* ************************************************************************** */
int8_t *GetStrArgv(int32_t argc,int8_t *argv[], int8_t *str)

{

 int32_t        idum;

 idum = FindOpt(argc,argv, str);

 if (idum != -1 && idum <= argc-2)
   {
    return argv[idum+1];
   }
 else
   {
    fprintf(stderr,"\n %s must be specified ! \n\n",str);
    exit(1);
   }
}
/* ************************************************************************** */
double GetDoubleArgv(int32_t argc,int8_t *argv[], int8_t *str)

{

 int32_t        idum;

 double     DoubleVal;

 idum = FindOpt(argc,argv, str);

 if (idum != -1 && idum <= argc-2)
   {
    sscanf(argv[idum + 1], "%lf", &DoubleVal);
    return DoubleVal;
   }
 else
   {
    fprintf(stderr,"\n %s must be specified ! \n\n",str);
    exit(1);
   }
}
int32_t GetArgEnvGraceVersion(int32_t argc, int8_t *argv[])
/*----------------------------------------------------------------------------->
/ purpose: return version number based on argument list with key -version
/          or if arg list is not defined then extract version from enviroment
/          variable GRACE_PRODUCT_VERSION
/
/ coded by: Gerhard L.H. Kruizinga                           06/25/2001
/
/ return: version number 0<= <=99 or -1L on failure
/
/ note:   version number defaults to 0 if neither argv or GRACE_PRODUCT_VERSION
/         is defined
/
<-----------------------------------------------------------------------------*/
{
 int32_t VersionNumber,arg_ndx;

 int8_t VersionChar[1000];

 VersionNumber = 0;               

 if (getenv("GRACE_PRODUCT_VERSION"))
 {
   arg_ndx = FindOpt(argc,argv, "-version");
   if (arg_ndx != -1)
   {
     VersionNumber = GetLongArgv(argc,argv,"-version");
   }
   else
   {
     strcpy(VersionChar,getenv("GRACE_PRODUCT_VERSION"));
     if (sscanf(VersionChar,"%d",&VersionNumber) != 1)
     {
       fprintf(stderr,
            "\n VersionNumber in enviroment variable GRACE_PRODUCT_VERSION = \"%s\"\n"
            ,VersionChar);
       fprintf(stderr," must be a numeral!! \n\n");
       return -1L;
     }   
     if (VersionNumber < 0 || VersionNumber > 99)
     {   
       fprintf(stderr,   
        "\n Version number in enviroment variable GRACE_PRODUCT_VERSION = %s\n",
               VersionChar);
       fprintf(stderr," needs to be defined within range >0 and <= 99 !!\n\n");
       return -1L;             
     }   
   }     
 }
 else    
 {                          
   arg_ndx = FindOpt(argc,argv, "-version");
   if (arg_ndx != -1) VersionNumber = GetLongArgv(argc,argv,"-version");
 }
 
 if (VersionNumber < 0 || VersionNumber > 99)
 {
   fprintf(stderr,"\n Specified Version number  = %d ", VersionNumber); 
   fprintf(stderr," needs to be defined within range >0 and <= 99 !!\n\n");
   exit(1);
 }

 return VersionNumber;
}
