/* $Id: GRACEio_prototypes.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */

#ifndef _GRACEio_prototypes_h_
#define _GRACEio_prototypes_h_

#include "GRACEdefs.h"
#include <string.h>

/*----------------------------------------------------------------------------->
/ Function Prototypes 
<-----------------------------------------------------------------------------*/

int32_t DecodeTMFileName(
         int8_t *TMFilename ,/*Standard Telemetry Filename                      */ 
         int8_t *Satellite,  /*Satellite indicator (A = GRACE A and B = GRACE B"*/
         int8_t *Version,    /*Version string of file for Date                  */
         int8_t *DateString, /*Timestamp in yyyy-mm-dd hh:mm:ss                 */
         int32_t *DateT2000,  /*Timestamp in seconds past 01/01/2000 12:00:00    */
         int8_t *Station,    /*Downlink Station Name (eg NZ = Neustrelitz)      */
         int8_t *ProcCenter, /*Processing center Name (eg RDC)                  */
         int8_t *DataStream, /*Satellite Data Stream (eg RT)                    */
         int8_t *Datatype    /*Date type (eg HK SC)                             */
         );

/*----------------------------------------------------------------------------*/
void ConstructFileName(
         int8_t *Satellite, /* Satellite indicator (A = GRACE A and B = GRACE B */
         double Time,     /* time for filename time tag (sec past 2000)       */
         int32_t VersionNumber, /* Version number of data product                */
         int8_t *FileTypeName, /* Standard file type name (aka filekey,eg ACC1B)*/
         int8_t *Filename,  /* output filename                                  */
         int8_t *ext        /* file extension                                   */
         );
/*----------------------------------------------------------------------------*/
void ConstructFileNameVersion(
         int8_t *Satellite, /* Satellite indicator (A = GRACE A and B = GRACE B */
         double Time,     /* time for filename time tag (sec past 2000)       */
         int8_t *Version,   /* Version string                                   */
         int8_t *FileTypeName, /* Standard file type name (aka filekey,eg ACC1B)*/
         int8_t *Filename,  /* output filename                                  */
         int8_t *ext        /* file extension                                   */
         );

void LinkTime 
               (
                 int8_t *linktimelabel     /* pointer to link time label        */    
               );

/*----------------------------------------------------------------------------*/
void GetBaseFilename 
               (int8_t *filename,          /* filename (including path)         */
                int8_t *BaseFilename       /* filename without path             */
               );

/*----------------------------------------------------------------------------*/
int32_t GetFileTypeName
               (
                int32_t SectorPointer,       /* file type pointer                */
                int8_t *filetype            /* Pointer to Filetype txt          */
               );

/*----------------------------------------------------------------------------*/
int32_t GetAcc1aProdName
               (
                int32_t SectorPointer,       /* Acc1a prodcuct pointer           */
                int8_t *filetype            /* Pointer to Acc1a product txt     */
               );

/*----------------------------------------------------------------------------*/

int32_t GetHeaderLabelName
               (
                int32_t HeaderRecPointer,    /* header record label pointer      */
                int8_t *filetype            /* Pointer to Filetype txt          */
               );

/*----------------------------------------------------------------------------*/

int32_t GetAcc1aProd
               (
                int8_t *filetype            /* Pointer to Acc1a product txt     */
               );

/*----------------------------------------------------------------------------*/

int32_t GetFileType
               (
                int8_t *filetype            /* Pointer to Filetype txt          */
               );

/*----------------------------------------------------------------------------*/

int32_t GetHeaderLabel
               (
                int8_t *filetype            /* Pointer to Filetype txt          */
               );

/*----------------------------------------------------------------------------*/
void InitializeHeaders();
/*----------------------------------------------------------------------------*/
void InitializeHeaderStruct(FileHeader_t *header);
/*----------------------------------------------------------------------------*/


double GetUTC2000TimeTag
               (
               );
void GetUTCTimeTag
               (
                int8_t            *time_tag /* Pointer to array with UTC timetag*/
               );
/*----------------------------------------------------------------------------*/
boolean WriteFileHeader 
	       (
/* input */
		FILE            *dst,     /* Pointer to data file             */
		FileHeader_t    *header   /* Pointer to header struct         */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/
/*----------------------------------------------------------------------------*/
boolean ReadFileHeader 
	       (
/* input */
		FILE            *src,     /* Pointer to data file             */
/* output */
		FileHeader_t    *header   /* Pointer to header struct         */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/
/*----------------------------------------------------------------------------*/
boolean PutFileHeader 
	       (
/* input */
                FILE *dst,                /* pointer to destination file      */
                FileHeader_t *header,     /* pointer to header struct         */
                int32_t FileType,            /* file type integer pointer        */
                int32_t FormatType,          /* format type binary/ascii         */
                int32_t NHeaderRecord,       /* number of header records         */
                int32_t init_flag,           /* initialization flag 0,1          */
                int8_t *ProdAgency,         /* producer agency label            */
                int8_t *ProdInstitution,    /* producer institution label       */
                int8_t *SoftwareVersion,    /* software version label           */
                int8_t *SoftwareLinkTime,   /* software version label           */
                int8_t *Documentation,      /* documentation label              */
                int8_t *SatName,            /* satellite name label             */
                int8_t *SensorName,         /* sensor name label                */
                int8_t *TimeEpoch,          /* Time epoch label                 */
                double TfirstObs,         /* Time first observation           */
                double TlastObs,          /* Time last observation            */
                int32_t   Nobs,              /* number of observations           */
                int8_t *TimeTag,            /* time tag label FIRST,FINAL       */
                int8_t *FileName,           /* Filename                         */
                int8_t *ProcessLevel        /* Processing Level (1A or 1B)      */
               );
/* return:	1	normal return
                0	End Of File reached 
*/
/*----------------------------------------------------------------------------*/
boolean ModifyFileHeader 
	       (
/* input */
		FileHeader_t    *header,   /* Pointer to header struct        */
                int32_t            RecPointer,/* Pointer to header record        */
                int8_t            *ModString /* Pointer to string to be inserted*/
                                           /* in header                       */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
boolean ReadGFD1XFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
/* output */
		GFD1X_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GFD1X_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteGFD1XFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
		GFD1X_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GFD1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintGFD1XFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		GFD1X_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GFD1X_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean WrAsciiGFD1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to flight data file      */
               GFD1X_t          *record   /* Pointer to flight data record    */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadGPS1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
/* output */
		GPS1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteGPS1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
		GPS1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean PrintGPS1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		GPS1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadGPS1BFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
/* output */
		GPS1B_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1B_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteGPS1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
		GPS1B_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1B_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean PrintGPS1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		GPS1B_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (GPS1B_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean WrAsciiGPS1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to flight data file      */
               GPS1B_t          *record   /* Pointer to flight data record    */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadKBR1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
/* output */
		KBR1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (KBR1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteKBR1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Flight Data Format*/
                                          /* file                             */
		KBR1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (KBR1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean PrintKBR1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		KBR1A_t       	*record   /* Pointer to GPS Flight Data Format*/
                                          /* struct (KBR1A_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadGNV1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Navigation 1A Data*/
                                          /* Format file                      */
/* output */
		GNV1A_t       	*record   /* Pointer to GPS Navigation 1A Data*/
                                          /* Format struct (GNV1A_t)          */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteGNV1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Navigation 1A Data*/
                                          /* Format file                      */
		GNV1A_t       	*record   /* Pointer to GPS Navigation 1A Data*/
                                          /* Format struct (GNV1A_t)          */
	       );

/*----------------------------------------------------------------------------*/
void PrintGNV1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		GNV1A_t       	*record   /* Pointer to GPS Navigation 1A Data*/
                                          /* Format struct (GNV1A_t)          */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiGNV1AFRecord
               (
/* input */
               FILE             *src,     /* Pointer to GPS Navigation 1A file*/
               GNV1A_t          *record   /* Pointer to GPS Nav. 1B record    */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadSCA1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to SCA 1A Data Format    */
                                          /* file                             */
/* output */
		SCA1A_t       	*record   /* Pointer to SCA 1A Data Format    */
                                          /* struct (SCA1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteSCA1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to SCA 1A Data Format    */
                                          /* file                             */
		SCA1A_t       	*record   /* Pointer to SCA 1A Data Format    */
                                          /* struct (SCA1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiSCA1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to SCA 1A Data Format    */
                                          /* file                             */
		SCA1A_t       	*record   /* Pointer to SCA 1A Data Format    */
                                          /* struct (SCA1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintSCA1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		SCA1A_t       	*record   /* Pointer to SCA 1A Data Format    */
                                          /* struct (SCA1A_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean ReadACC1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to ACC 1A Data Format    */
                                          /* file                             */
/* output */
		ACC1A_t       	*record   /* Pointer to ACC 1A Data Format    */
                                          /* struct (ACC1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteACC1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to ACC 1A Data Format    */
                                          /* file                             */
		ACC1A_t       	*record   /* Pointer to ACC 1A Data Format    */
                                          /* struct (ACC1A_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean WriteAHK1XFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to AHK 1X Data Format    */
                                          /* file                             */
		AHK1X_t       	*record   /* Pointer to AHK 1X Data Format    */
                                          /* struct (AHK1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintACC1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		ACC1A_t       	*record   /* Pointer to ACC 1A Data Format    */
                                          /* struct (ACC1A_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean ReadGNV1BFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to GPS Navigation 1B Data*/
                                          /* Format file                      */
/* output */
		GNV1B_t       	*record   /* Pointer to GPS Navigation 1B Data*/
                                          /* Format struct (GNV1B_t)          */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteGNV1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to GPS Navigation 1B Data*/
                                          /* Format file                      */
		GNV1B_t       	*record   /* Pointer to GPS Navigation 1B Data*/
                                          /* Format struct (GNV1B_t)          */
	       );

/*----------------------------------------------------------------------------*/
void PrintGNV1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		GNV1B_t       	*record   /* Pointer to GPS Navigation 1B Data*/
                                          /* Format struct (GNV1B_t)          */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiGNV1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to GPS Navigation 1B file*/
               GNV1B_t          *record   /* Pointer to GPS Nav. 1B record    */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadTIM1XFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to TIM1X Data Format    */
                                          /* file                             */
/* output */
		TIM1X_t       	*record   /* Pointer to TIM1X Data Format    */
                                          /* struct (TIM1X_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/
/*----------------------------------------------------------------------------*/
boolean ReadPCI1AFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to PCI 1A Data Format    */
                                          /* file                             */
/* output */
		PCI1A_t       	*record   /* Pointer to PCI 1A Data Format    */
                                          /* struct (PCI1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WritePCI1AFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to PCI 1A Data Format    */
                                          /* file                             */
		PCI1A_t       	*record   /* Pointer to PCI 1A Data Format    */
                                          /* struct (PCI1A_t)                 */
	       );
/*----------------------------------------------------------------------------*/
boolean ReadSCA1BFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to SCA 1B Data Format    */
                                          /* file                             */
/* output */
		SCA1B_t       	*record   /* Pointer to SCA 1B Data Format    */
                                          /* struct (SCA1B_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteSCA1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to SCA 1B Data Format    */
                                          /* file                             */
		SCA1B_t       	*record   /* Pointer to SCA 1B Data Format    */
                                          /* struct (SCA1B_t)                 */
	       );
/*----------------------------------------------------------------------------*/
boolean WriteTIM1XFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to TIM1X  Data Format    */
                                          /* file                             */
		TIM1X_t       	*record   /* Pointer to TIM1X  Data Format    */
                                          /* struct (TIM1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintSCA1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		SCA1B_t       	*record   /* Pointer to SCA 1B Data Format    */
                                          /* struct (SCA1B_t)                 */
	       );
/*----------------------------------------------------------------------------*/
void PrintTIM1XFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		TIM1X_t       	*record   /* Pointer to TIM1X Data Format     */
                                          /* struct (TIM1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiPCI1AFRecord
               (
/* input */
               FILE             *src,     /* Pointer to star camera file      */
               PCI1A_t          *record   /* Pointer to star camera record    */
	       );
/*----------------------------------------------------------------------------*/
boolean WrAsciiSCA1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to star camera file      */
               SCA1B_t          *record   /* Pointer to star camera record    */
	       );


/*----------------------------------------------------------------------------*/
boolean WrAsciiTIM1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to OBDH mapping file     */
               TIM1X_t          *record   /* Pointer to OBDH mapping record   */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadKBR1BFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to KBR 1B Data Format    */
                                          /* file                             */
/* output */
		KBR1B_t       	*record   /* Pointer to KBR 1B Data Format    */
                                          /* struct (KBR1B_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteKBR1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to KBR 1B Data Format    */
                                          /* file                             */
		KBR1B_t       	*record   /* Pointer to KBR 1B Data Format    */
                                          /* struct (KBR1B_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintKBR1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		KBR1B_t       	*record   /* Pointer to KBR 1B Data Format    */
                                          /* struct (KBR1B_t)                 */
	       );
/*----------------------------------------------------------------------------*/
boolean WrAsciiKBR1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to KBR Level 1B file     */
               KBR1B_t          *record   /* Pointer to KBR Level 1B x record */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadACC1BFRecord 
	       (
/* input */
		FILE		*src,     /* Pointer to ACC 1B Data Format    */
                                          /* file                             */
/* output */
		ACC1B_t       	*record   /* Pointer to ACC 1B Data Format    */
                                          /* struct (ACC1B_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteACC1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to ACC 1B Data Format    */
                                          /* file                             */
		ACC1B_t       	*record   /* Pointer to ACC 1B Data Format    */
                                          /* struct (ACC1B_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintACC1BFRecord 
	       (
/* input */
		FILE		*dst,     /* Pointer to output file           */
		ACC1B_t       	*record   /* Pointer to ACC 1B Data Format    */
                                          /* struct (ACC1B_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean WrAsciiACC1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to Accelerometer file    */
               ACC1B_t          *record   /* Pointer to Accelerometer record  */
	       );


/*----------------------------------------------------------------------------*/
boolean ReadXXXVOFRecord 
	       (
/* input */
		FILE	          *src,   /* Pointer to XXXVOB Data Format    */
                                          /* file                             */
/* output */
		XXXVO_t           *record /* Pointer to XXXVOB Data Format    */
                                          /* struct (XXXVO_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteXXXVOFRecord 
	       (
/* input */
		FILE		  *dst,   /* Pointer to XXXVOB Data Format    */
                                          /* file                             */
		XXXVO_t           *record /* Pointer to XXXVOB Data Format    */
                                          /* struct (XXXVO_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintXXXVOFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file          */
		XXXVO_t          *record  /* Pointer to XXXVOB Data Format   */
                                          /* struct (XXXVO_t)                */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiXXXVOFRecord
               (
/* input */
               FILE             *src,     /* Pointer to vector orientation    */
                                          /* file                             */
               XXXVO_t          *record   /* Pointer to vector orientation    */
                                          /* record                           */
	       );

/*----------------------------------------------------------------------------*/

boolean ReadIOA1BFRecord 
	       (
/* input */
		FILE	         *src,    /* Pointer to IOA 1B Data Format    */
                                          /* file                             */
/* output */
		IOA1B_t          *record  /* Pointer to IOA 1B Data Format    */
                                          /* struct (IOA1B_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
boolean WriteIOA1BFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to IOA 1B Data Format    */
                                          /* file                             */
		IOA1B_t          *record  /* Pointer to IOA 1B Data Format    */
                                          /* struct (IOA1B_t)                 */
	       );

/*----------------------------------------------------------------------------*/
void PrintIOA1BFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		IOA1B_t          *record  /* Pointer to IOA 1B Data Format    */
                                          /* struct (IOA1B_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean WrAsciiIOA1BFRecord
               (
/* input */
               FILE             *src,     /* Pointer to accelerometer inertial*/
                                          /* orientation file                 */
               IOA1B_t          *record   /* Pointer to accelerometer inertial*/
                                          /* orientation record               */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadOSCFQFRecord 
	       (
/* input */
		FILE	         *src,    /* Pointer to OSCFQ Data Format     */
                                          /* file                             */
/* output */
		OSCFQ_t          *record  /* Pointer to OSCFQ Data Format     */
                                          /* struct (OSCFQ_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/

/*----------------------------------------------------------------------------*/
void PrintOSCFQFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		OSCFQ_t          *record  /* Pointer to OSCFQ Data Format     */
                                          /* struct (OSCFQ_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean WrAsciiOSCFQFRecord
               (
/* input */
               FILE             *src,     /* Pointer to USO file              */
               OSCFQ_t          *record   /* Pointer to USO record            */
	       );

/*----------------------------------------------------------------------------*/
boolean WriteOSCFQFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to OSCFQ  Data Format    */
                                          /* file                             */
		OSCFQ_t          *record  /* Pointer to OSCFQ  Data Format    */
                                          /* struct (OSCFQ_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadMAG1XFRecord 
	       (
/* input */
		FILE		 *src,    /* Pointer to MAG Data Format      */
                                          /* file                            */
		MAG1X_t          *record  /* Pointer to MAG Data Format      */
                                          /* struct (MAG1X_t)                */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintMAG1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file          */
		MAG1X_t          *record  /* Pointer to MAG Data Format      */
                                          /* struct (MAG1X_t)                */
	       );



/*----------------------------------------------------------------------------*/
boolean WrAsciiMAG1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to MAG file             */
               MAG1X_t          *record   /* Pointer to MAG record           */
	       );


/*----------------------------------------------------------------------------*/
boolean WriteMAG1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to MAG Data Format       */
                                          /* file                             */
		MAG1X_t          *record  /* Pointer to MAG Data Format       */
                                          /* struct (MAG1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadIHK1XFRecord 
	       (
/* input */
		FILE		 *src,    /* Pointer to IPU HK   Data Format  */
                                          /* file                             */
		IHK1X_t          *record  /* Pointer to IPU HK   Data Format  */
                                          /* struct (IHK1X_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintIHK1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		IHK1X_t          *record  /* Pointer to IPU HK   Data Format  */
                                          /* struct (IHK1X_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean WrAsciiIHK1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to IPU HK   file         */
               IHK1X_t          *record   /* Pointer to IPU HK   record       */
	       );

/*----------------------------------------------------------------------------*/
boolean WriteIHK1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to IPU HK   Data Format  */
                                          /* file                             */
		IHK1X_t          *record  /* Pointer to IPU HK   Data Format  */
                                          /* struct (THR1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadTHR1XFRecord 
	       (
/* input */
		FILE		 *src,    /* Pointer to Thruster Data Format  */
                                          /* file                             */
		THR1X_t          *record  /* Pointer to Thruster Data Format  */
                                          /* struct (THR1X_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintTHR1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		THR1X_t          *record  /* Pointer to Thruster Data Format  */
                                          /* struct (THR1X_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean WrAsciiTHR1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to Thruster file         */
               THR1X_t          *record   /* Pointer to Thruster record       */
	       );

/*----------------------------------------------------------------------------*/
boolean WriteTHR1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to Thruster Data Format  */
                                          /* file                             */
		THR1X_t          *record  /* Pointer to Thruster Data Format  */
                                          /* struct (THR1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadTNK1XFRecord      
               (
/* input */
                FILE             *src,    /* Pointer to Tank Data Format      */
                                          /* file                             */
                TNK1X_t          *record  /* Pointer to Tank Data Format      */
                                          /* struct (TNK1X_t)                 */
               );
/* return:      1       normal return
                0       End Of File reached
*/             
                 
               
/*----------------------------------------------------------------------------*/
void PrintTNK1XFRecord
               (    
/* input */    
                FILE             *dst,    /* Pointer to output file           */
                TNK1X_t          *record  /* Pointer to Tank Data Format      */
                                          /* struct (TNK1X_t)                 */
               );

               
                
/*----------------------------------------------------------------------------*/
boolean WrAsciiTNK1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to Tank file             */
               TNK1X_t          *record   /* Pointer to Tank Data record      */
               );
                
/*----------------------------------------------------------------------------*/
boolean WriteTNK1XFRecord
               (
/* input */
                FILE             *dst,    /* Pointer to Tank Data Format      */
                                          /* file                             */
                TNK1X_t          *record  /* Pointer to Tank Data Format      */
                                          /* struct (TNK1X_t)                 */
               );


/*----------------------------------------------------------------------------*/
boolean ReadMAS1XFRecord 
	       (
/* input */
		FILE		 *src,    /* Pointer to Mass Data Format      */
                                          /* file                             */
		MAS1X_t          *record  /* Pointer to Mass Data Format      */
                                          /* struct (MAS1X_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintMAS1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		MAS1X_t          *record  /* Pointer to Mass Data Format      */
                                          /* struct (MAS1X_t)                 */
	       );



/*----------------------------------------------------------------------------*/
boolean WrAsciiMAS1XFRecord
               (
/* input */
               FILE             *src,     /* Pointer to Mass file             */
               MAS1X_t          *record   /* Pointer to Mass Data record      */
	       );

/*----------------------------------------------------------------------------*/
boolean WriteMAS1XFRecord 
	       (
/* input */
		FILE		 *dst,    /* Pointer to Mass Data Format      */
                                          /* file                             */
		MAS1X_t          *record  /* Pointer to Mass Data Format      */
                                          /* struct (MAS1X_t)                 */
	       );

/*----------------------------------------------------------------------------*/
boolean ReadHRT1XFRecord(
/* input */
               FILE *src,                /* Pointer to HRT file               */ 
               HRT1X_t *record           /* Pointer to HRT struct file        */
               );
/* return:	1	normal return
                0	End Of File reached 
*/

boolean WrAsciiHRT1XFRecord(FILE *dst, HRT1X_t *record);
boolean WriteHRT1XFRecord(FILE *dst, HRT1X_t *record);
void  PrintHRT1XFRecord(FILE *dst, HRT1X_t *record);

/*----------------------------------------------------------------------------*/
boolean ReadAHK1AFRecord  
/* Not yet implemented */
	       (
/* input */
		FILE		 *src,    /* Pointer to Accel. HK Data Format */
                                          /* file                             */
		AHK1X_t          *record  /* Pointer to Accel. HK Data Format */
                                          /* struct (AHK1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintAHK1AFRecord 
/* Not yet implemented */
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		AHK1X_t          *record  /* Pointer to Accel. HK data Format */
                                          /* struct (AHK1A_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean WriteAHK1AFRecord 
/* Not yet implemented */
	       (
/* input */
		FILE		 *dst,    /* Pointer to Accel. HK Data Format */
                                          /* file                             */
		AHK1X_t          *record  /* Pointer to Accel. HK Data Format */
                                          /* struct (AHK1A_t)                 */
	       );
/*----------------------------------------------------------------------------*/
boolean ReadIHK1AFRecord 
/* Not yet implemented */
	       (
/* input */
		FILE		 *src,    /* Pointer to IPU HK Data Format    */
                                          /* file                             */
		IHK1X_t          *record  /* Pointer to IPU HK Data Format    */
                                          /* struct (IHK1A_t)                 */
	       );
/* return:	1	normal return
                0	End Of File reached 
*/


/*----------------------------------------------------------------------------*/
void PrintIHK1AFRecord 
/* Not yet implemented */
	       (
/* input */
		FILE		 *dst,    /* Pointer to output file           */
		IHK1X_t          *record  /* Pointer to IPU HK data Format    */
                                          /* struct (IHK1A_t)                 */
	       );


/*----------------------------------------------------------------------------*/
boolean WriteIHK1AFRecord 
/* Not yet implemented */
	       (
/* input */
		FILE		 *dst,    /* Pointer to IPU HK Data Format    */
                                          /* file                             */
		IHK1X_t          *record  /* Pointer to IPU HK Data Format    */
                                          /* struct (IHK1A_t)                 */
	       );


/*----------------------------------------------------------------------------*/
int32_t VerifyFilename (

/* input */
                int8_t *filename,         /* input filename to be verified      */
                int8_t *BaseFilename,     /* input filename with path removed   */
                int8_t *FileKey,          /* retrieved filekey (valid file only)*/
                int32_t *year,              /* retrieved year (valid file only)   */
                int32_t *month,             /* retrieved month (valid file only)  */
                int32_t *day,               /* retrieved day  (valid file only)   */
                int8_t *SatId,            /* retrieved satid (valid file only)  */
                int32_t *VersionNumber,    /* retrieved vs No. (valid file only) */
                int8_t *VersionChar,      /* retrieved vs Str (valid file only) */
                int8_t *ext               /* retrieved ext (valid file only)    */
                );
/*----------------------------------------------------------------------------*/
void MakeL1BOutputFilename (
                int8_t *l1a_inputfile,      /* input l1a filename               */
                int8_t *l1b_inputfilename,  /* input l1b filename or null       */
                int8_t *l1b_outputfilename, /* output l1b filename              */
                int8_t *l1b_reportfilename, /* output l1b report filename       */
                int8_t *l1b_charversion,    /* version string for l1b filenames */
                                          /* in case of standard filename     */
                int8_t *l1b_filesat         /* satellite id, ignore if set to \0*/
                );
/*----------------------------------------------------------------------------*/
int32_t little_endian();
/*----------------------------------------------------------------------------*/
int32_t swapbyte(
                int8_t buf[],               /* byte array to be swapped         */
                size_t num_bytes          /* size of variable to be swapped   */
             );
/*----------------------------------------------------------------------------*/
size_t fread_grace(
                void *ptr, 
                size_t  size,  
                size_t  nitems,  
                FILE *stream
             );  
/*----------------------------------------------------------------------------*/
size_t fwrite_grace(
                void *ptr, 
                size_t  size,  
                size_t  nitems,  
                FILE *stream
             );  
boolean ReformatRCStag(int8_t *SoftwareVersion,    /* RCS software version label*/
                       int8_t *NewLabel            /* Processing level (1A or 1B)*/
                      );
boolean LoadClockFile(FILE *clk,double **xt, double **yt,int32_t *Np);
void AddTSsuppid2qualflg(signed char   TSsuppId, uint8_t *qualflg);
int32_t Check4GPSinQualflag(uint8_t qualflg);
boolean ReadILG1XFRecord(FILE *src, ILG1X_t *record);
boolean WriteILG1XFRecord(FILE *dst, ILG1X_t *record);
int32_t OBDH2RCVRtime(int32_t obdh_intg, int32_t obdh_frac, int32_t *gps_intg, int32_t *gps_frac,
                  int32_t icu_blk_nr, FILE *tim);
int32_t Write_CMT_command(CMT_command_t *record);
int32_t Write_SCA2K_command(SCA2K_command_t *record);
boolean AddInputFile2Header(FileHeader_t *header,int8_t *filekey, 
                            int8_t *input_filename,int8_t *input_file_ttag,
                            int8_t *version, int8_t *linktime);
boolean CopyInputFile2Header(FileHeader_t *src_header,int8_t *filekey, 
                             FileHeader_t *dst_header);




#endif	/* _GRACEio_prototypes_h_ */
