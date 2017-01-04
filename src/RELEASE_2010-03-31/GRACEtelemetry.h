/* $Id: GRACEtelemetry.h,v 1.1 2009/06/03 22:53:17 glk Exp glk $ */
#include <string.h>
#include "GRACEfiletype.h"

#ifndef _GRACEtelemetry_h_ 
#define _GRACEtelemetry_h_ 

#define MAXPACKET           500
#define MAXNRPACKETS       (MAXPACKET+1) 
#define PKT_NOTDEFINED       -1
#define NOPKTPOINTER         -1
#define PACKETCHARMAX       100
#define NRTMSTATS             4
#define IPTOTALPKT            0
#define IPGOODPKT             1
#define IPBADPKT              2
#define IPBADBYTES            3
#define IPINVALIDPACKET       MAXPACKET
#define NRICUCOEF             4
#define NRICUCALTEMPS         3
#define NRICUBOARDS           2
#define NRICUSATS             2
#define IPICUBOARDNOM         0
#define IPICUBOARDRED         1
#define IPICUGRACEA           0
#define IPICUGRACEB           1

#define IPPREVTS              0
#define IPCURRTS              1
#define IPNEXTTS              2
#define MAXICUPKTS            100

#define MAXHKCALCHAR          200
#define MAXHKCALLINES         200
#define MAXPOLYNORDER         10
#define NOT_USED              0L
#define USE                   1L
#define MAXHKSAT              2L

#define NTHRUST_MAX           100
#define MAXTHRUSTERS          14
#define IPGRACEA              0
#define IPGRACEB              1

#define NSOESENSORS          28
#define IPSOEIPU              0
#define IPSOEKBR              1
#define IPSOESCA              2
#define IPSOEACC              3
#define IPSOEUSO              4
#define IPSOEVGN              5
#define IPSOEVGB              6
#define IPSOEVGO              7
#define IPSOEVSL              8
#define IPSOEVCM              9
#define IPSOEVKB             10
#define IPSOEQSA             11
#define IPSOEQSB             12
#define IPSOEMTE1            13
#define IPSOEMTE2            14
#define IPSOEQKS             15
#define IPSOEIPUR            16
#define IPSOEK_MI            17
#define IPSOEKAMI            18
#define IPSOEKTOFF           19
#define IPSOEAOCS            20
#define IPSOEICUVP           21
#define IPSOEACCR            22
#define IPSOEACCT            23
#define IPSOEMANV            24
#define IPSOECMCAL           25
#define IPSOEKBRCAL          26
#define IPSOEOCC             27

#define GPSTIMEOFFSET -630763200

#define MAXAPPDATA 65535 /* this number is determined by adding all max packet*/
                         /* lengths in file GRACE_TM_properties.txt. This is  */
                         /* the maximum amount of bytes that can occur before */
                         /* a time stamp packet must be found in the byte     */
                         /* stream. 65535 is chose because max uint16_t */

#define NTHRUSTER  7
#define MAXPARAM  10
#define MAXWORD  100

typedef struct TMpackets_t            /* struct for TM packet properties      */
        {   
          int8_t name[PACKETCHARMAX];   /* TM packet name                       */
          int32_t PacketId;              /* Packet Id                            */
          int32_t MinSuppPacketId;       /* Minimum Supplemental PacketId        */
          int32_t MaxSuppPacketId;       /* Minimum Supplemental PacketId        */
          int32_t MinPacketLength;       /* Minimum Packet length (bytes)        */
          int32_t MaxPacketLength;       /* Maximum Packet Length (bytes)        */
        } TMpackets_t;

typedef struct TimeStamp_t            /* struct for data in Time Stamp packet */
        {
          uint8_t *PacketId;     /* Timestamp packet Id                  */
          signed char   *SuppPacketId;   /* Timestamp supplementary packet Id    */
          uint16_t *Dlen;        /* Timestamp data packet length         */
          int32_t *Time;                  /* time tag (GPS sec past 01/06/80 0:0:0*/
        } TimeStamp_t;

typedef struct ApplicationPkt_t        /* struct for data in any application   */
                                       /* packet                               */
        {
          uint8_t *PacketId;     /* Application packet Id                */
          signed char   *SuppPacketId;   /* Application supplementary packet Id  */
          uint16_t *Dlen;        /* Application data packet length       */
          uint8_t *Data;         /* pointer to data segment in App packet*/
        } ApplicationPkt_t;

typedef struct ICUDAcoef_t             /* coefficient for ICU digital to analog*/
                                       /* data conversion                      */
        {
          int8_t sccsid[MAXHKCALCHAR];   /* sccs id of set up file               */
          double temps[NRICUCALTEMPS]; /* temperature for coefficients in deg C*/
          double coeff[NRACC1APODS][NRICUCOEF][NRICUCALTEMPS]; /* coefficients */
        } ICUDAcoef_t;

typedef struct ICU_packet_t            /* struct containing ICU packet info    */
        {
          uint16_t pkt_length;   /* icu packet length (including header  */
          int8_t     pkt_name;           /* icu packet name                      */
          int8_t     pkt_chksum_status;  /* icu checksum status                  */
                                       /* 0 = OK, 1 = invalid                  */
          uint16_t block_nr;     /* icu packet block number              */ 
        } ICU_packet_t;

typedef struct ICUAP_info_t            /* struct containing ICU application    */
                                       /* packet information                   */
        {
          uint8_t AP_pktId;      /* Application Packet Id                */
          int8_t AP_SuppId;              /* Application Packet Supplementary Id  */
          uint16_t AP_Dlen;      /* Application Packet Length            */
          int32_t TimeStamp[3];           /* previous, current and next timestamp */
          int8_t TS_SuppId[3];           /* prev, current and next timestamp     */
                                       /* supplementary ID where               */
                                       /* 0 = SCET (space craft elapsed time)  */
                                       /* 1 = GPS time                         */
                                       /* 2 = SCET + PulseSync                 */
                                       /* 3 = GPS time + PulseSync             */
          int32_t NicuPkts;               /* number of ICU pkts in ICU AP packet  */
          ICU_packet_t pkt[MAXICUPKTS];/* icu packet info for all ICU packets  */
          uint16_t blocknrs[3];  /* previous and next ICU block numbers  */
        } ICUAP_info_t;

typedef struct TM_MAG_t     /*            Magnetometer HK packet struct        */
        {
          int16_t Updated;    /*1 indicates record data updated                  */
          int16_t MfvX_RAW;   /*x-axis component of measured earth magnetic field*/
          int16_t MfvY_RAW;   /*y-axis component of measured earth magnetic field*/
          int16_t MfvZ_RAW;   /*z-axis component of measured earth magnetic field*/
          int16_t Torque1A;   /*current of magnetorquer 1 A (positive current x) */
          int16_t Torque2A;   /*current of magnetorquer 2 A (positive current y) */
          int16_t Torque3A;   /*current of magnetorquer 3 A (positive current z) */
          int16_t Torque1B;   /*current of magnetorquer 1 B (negative current x) */
          int16_t Torque2B;   /*current of magnetorquer 2 B (negative current y) */
          int16_t Torque3B;   /*current of magnetorquer 3 B (negative current z) */
          int16_t MF_BCalX;   /*mag field calibration factors (x)                */
          int16_t MF_BCalY;   /*mag field calibration factors (y)                */
          int16_t MF_BCalZ;   /*mag field calibration factors (z)                */
          int16_t MT_Kr;      /*mag torquer calibration factor                   */
        } TM_MAG_t;
typedef struct TM_ANARAW_SCT_t /* Standard Calibration Thermistor(SCT) packet  */
        {
          int16_t value00;            /* Reference value for next 3 values       */ 
          int16_t value01;            /* Harness to MWA electronics +y red.      */
          int16_t value02;            /* ACC thermal Cage +z red                 */
          int16_t value03;            /* RFEA                                    */
          int16_t value04;            /* Reference value for next 15 values      */
          int16_t value05;            /* I/F Support Structure to MEP -y,red.    */
          int16_t value06;            /*ACC thermal Cage -z red.                 */
          int16_t value07;            /*OBDH                                     */
          int16_t tank_adap_negx;     /*Tank  adapter -x                         */
          int16_t value09;            /*Thruster Valve x ATH 12, 14, OTH 11, 21  */
          int16_t value10;            /*Platform nadir +y front                  */
          int16_t value11;            /*Platform nadir -y mid                    */
          int16_t value12;            /*Thruster Valve front y                   */
          int16_t value13;            /*Thruster Valve front z                   */
          int16_t tank_skin_negx;     /*Tank sphere -x                           */
          int16_t value15;            /*Solar Array -y side +x                   */
          int16_t value16;            /*MLI Heater-y                             */
          int16_t value17;            /*MLI Heater +y                            */
          int16_t value18;            /*MLI Heater MEP                           */
          int16_t value19;            /*Harness to MWA electronics -y red.       */
          int16_t value20;            /*Reference value for next 15 analog values*/
          int16_t value21;            /*I/F Support Structure to MEP +y, red.    */
          int16_t value22;            /*Harness to ACC Sensor red.               */
          int16_t value23;            /*Battery +/-y                             */
          int16_t value24;            /*PCDU                                     */
          int16_t value25;            /*Thruster Valve -y ATH 15,16,11,23        */
          int16_t value26;            /*Thruster Valve +x ATH 24, 22             */
          int16_t value27;            /*Platform nadir y front                   */
          int16_t value28;            /*I/F GPS Zenith Antenna                   */
          int16_t value29;            /*Thruster Valve front +y                  */
          int16_t value30;            /*Thruster Valve front +z                  */
          int16_t tank_skin_posx_red; /*Tank sphere +x redundant                 */
          int16_t value32;            /*Solar Array +y side +x                   */
          int16_t value33;            /*MLI Heater-y, redundant                  */
          int16_t value34;            /*MLI Heater +y, redundant                 */
          int16_t value35;            /*MLI Heater MEP, redundant                */
          int16_t value36;            /*Reference value of next 12 values        */
          int16_t value37;            /*RF-Sampling Unit red.                    */
          int16_t value38;            /*Battery +/-y redundant                   */
          int16_t tank_adap_posx;            /*Tank adapter +x                   */
          int16_t value40;            /*Thruster Valve +y ATH 26, 25 13, 21      */
          int16_t value41;            /*Magnetometer on Boom                     */
          int16_t value42;            /*Platform nadir +y mid                    */
          int16_t value43;            /*Spare                                    */
          int16_t value44;            /*Spare  candidate: IMU TBD,               */
          int16_t tank_skin_posx;     /*Tank sphere +x                           */
          int16_t tank_skin_negx_red; /*Tank sphere -x redundant                 */
          int16_t value47;            /*+x Panel Inside                          */
          int16_t value48;            /*I/F Support Structure to MEP mid red     */ 
          int16_t SctR0;              /*calibration parameter R0                 */
          int16_t SctRref;            /*calibration parameter Rref               */
        } TM_ANARAW_SCT_t;

typedef struct TM_PWR_PCDUHK_t      /* PCDU HK data packet                     */
        {
          int16_t value00;            /* Reference Voltage                       */
          int16_t value01;            /* (First) Solar Array Input Current 1 Main*/
          int16_t value02;            /* (First) Solar Array Input Current 1 Red */
          int16_t value03;            /* OBDH 1 Current                          */
          int16_t value04;            /* OBDH 2 Current                          */
          int16_t value05;            /* Receiver 1 Current                      */
          int16_t value06;            /* Receiver 2 Current                      */
          int16_t value07;            /* Transmitter Main Current                */
          int16_t value08;            /* Transmitter Red. Current                */
          int16_t value09;            /* MTE Main Current                        */
          int16_t value10;            /* MTE Red. Current                        */
          int16_t value11;            /* USO Right Current                       */
          int16_t value12;            /* USO Left Current                        */
          int16_t value13;            /* PARSU IPU Main Current                  */
          int16_t value14;            /* PARSU IPU Red. Current                  */
          int16_t value15;            /* SuperSTAR Accelerometer Current ICU Main*/
          int16_t value16;            /* SuperSTAR Accelerometer Current ICU Red.*/
          int16_t value17;            /* Getter Pump Current                     */
          int16_t value18;            /* PARSU MWA1 Current                      */
          int16_t value19;            /* (First Module) Heaters combined Current1*/
          int16_t value20;            /* Battery Charge Current Main             */
          int16_t value21;            /* Battery Charge Current Red.             */
          int16_t value22;            /* Battery Discharge Current Main          */
          int16_t value23;            /* Battery Discharge Current Red.          */
          int16_t value24;            /* Battery Temperature Main                */
          int16_t value25;            /* Battery Temperature Red.                */
          int16_t value26;            /* Battery Cell Pressure Main              */
          int16_t value27;            /* Battery Cell Pressure Red.              */
          int16_t value28;            /* Battery (Full) Voltage Main             */
          int16_t value29;            /* Battery (Full) Voltage Red.             */
          int16_t value30;            /* Battery 1 / 2 Voltage                   */
          int16_t value31;            /* PARSU MWA2 Current                      */
          int16_t regpress_1;         /* (First) Pressure Transducer Low 1       */
          int16_t tnkpress_1;         /* (First) Pressure Transducer High 1      */
          int16_t value34;            /* OC Thruster Combined Current            */
          int16_t value35;            /* (Second Module)Heaters combined Current2*/
          int16_t value36;            /* (Second)Solar Array Input Current 2 Main*/
          int16_t value37;            /* (Second)Solar Array Input Current 2 Red.*/
          int16_t regpress_2;         /* (Second) Pressure Transducer Low 2      */
          int16_t tnkpress_2;         /* (Second) Pressure Transducer High 2     */
          int16_t value40;            /* IMU (Gyro) Supply Main                  */
          int16_t value41;            /* AC Thruster Combined Current            */
          int16_t value42;            /* Analog HK Spare                         */
        } TM_PWR_PCDUHK_t;
typedef struct TM_AOC_STAT_t           /* AOCS status HK packet                */
        {
          int16_t NumErrs;               /* number of errors occured (wrap count)*/
          int16_t NumCmdsReceived;       /* number of cmds rcvd (wrap count)     */
          int16_t NumCmdsRejected;       /* number of cmds rejected (wrap count) */
          int16_t ActCtlFlags;           /* actuator control related flags       */
                                       /* b15: 1=>Thruster control by Ground;  */
                                       /*      0=>by OBSW                      */
                                       /* b14: 1=>Mag.Torquer control by Ground*/
                                       /*      0=>by OBSW                      */
                                       /* b13-b12: reserved                    */
                                       /* b11-b0: Thruster Enable bits         */
                                       /* b11 corresponds-- to THRUSTER_12,.., */
                                       /* b0 to THRUSTER_1                     */
          int32_t AccumThrOnTimes[12];    /* accum. thruster on-times [ms] (wrap) */ 
                                       /* for each 12 aoc thrusters            */
          int32_t AccumThrActivations[12];/* activation counters (wrap count)     */
                                       /* for each 12 aoc thrusters            */
          int32_t ActivationTimeStamp;    /* TimeStamp [s] of most recent thruster*/
                                       /* activation                           */
          int32_t AocModeChanges;         /* accum. number of mode changes (wrap) */
        } TM_AOC_STAT_t;

typedef struct tm_sc_stat_t            /* Timestamp + SC packet statistics     */
        {
          int32_t count;                  /* timestamp count in TM file           */
          int32_t TStime;                 /* timestamp time                       */
          int32_t TSstate;                /* timestamp suppid state               */
                                       /* 0 = SCET                             */
                                       /* 1 = GPS Time                         */
                                       /* 2 = SCET + Pulse sync                */
                                       /* 3 = GPS Time + Pulse sync            */
          int32_t TSstate_prev;           /* previous timestamp suppid state      */
                                       /* 0 = SCET                             */
                                       /* 1 = GPS Time                         */
                                       /* 2 = SCET + Pulse sync                */
                                       /* 3 = GPS Time + Pulse sync            */
          int32_t TSdt;                   /* Delta Time between current and       */
                                       /* previous TStime                      */
          int32_t IPUpktid;               /* IPU application packet Id            */
          int32_t IPUpktsize;             /* IPU application packet size          */
          int32_t ICUpktid;               /* ICU application packet Id            */
          int32_t ICUpktsize;             /* ICU application packet size          */
          int32_t NICUblocks;             /* Number of ICU datablock in ICU ap pkt*/
          int32_t DeltaBlockNr;           /* Delta Block number with first data   */
                                       /* block number of previous ICU ap pkt  */
          int32_t blocknrs[100];          /* ICU block numbers                    */
          int32_t blockgdel_flag[100];    /* ICU blocks for which gdel is set     */
                                       /* 0 = gdel flag not set, 1 = flag set  */
         } tm_sc_stat_t;

typedef struct hkcal_raw_t                /* struct containing raw hkcal info  */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_raw_t;

typedef struct hkcal_poly_t               /* struct containing poly hk cal info*/
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int32_t norder;                   /* order of polynomial               */
           double coeff[MAXPOLYNORDER];   /* polynomial coefficients           */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_poly_t;

typedef struct hkcal_scale_t               /* struct containing scale hkcal info*/
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           double scale;                  /* scale factor                      */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_scale_t;

typedef struct hkcal_f2args_t             /* struct containing f2 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f2 function         */
           int32_t arg3;                     /* argument 3 to f2 function (int32_t)  */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f2args_t;

typedef struct hkcal_f3args_t             /* struct containing f3 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f3 function         */
           int32_t arg3;                     /* argument 3 to f3 function (int32_t)  */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f3args_t;

typedef struct hkcal_f4args_t             /* struct containing f4 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg1[MAXHKCALCHAR];       /* argument 1 to f4 function         */
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f4 function         */
           int8_t arg3[MAXHKCALCHAR];       /* argument 2 to f4 function         */
           double arg4;                   /* argument 4 to f4 function (double)*/
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f4args_t;

typedef struct hkcal_f5args_t             /* struct containing f5 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg1[MAXHKCALCHAR];       /* argument 1 to f5 function         */
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f5 function         */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f5args_t;

typedef struct hkcal_f6args_t             /* struct containing f6 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg1[MAXHKCALCHAR];       /* argument 1 to f6 function         */
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f6 function         */
           int8_t arg3[MAXHKCALCHAR];       /* argument 2 to f6 function         */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f6args_t;

typedef struct hkcal_f7args_t             /* struct containing f7 arguments    */
         {
           int8_t on;                       /* function is active on = 1 else = 0*/
           int8_t arg1[MAXHKCALCHAR];       /* argument 1 to f7 function         */
           int8_t arg2[MAXHKCALCHAR];       /* argument 2 to f7 function         */
           int32_t npar;                     /* # of parameters need for function */
           int32_t par_ndx[MAXPOLYNORDER];   /* look up indeces for parameters    */
                                          /* used in function                  */
         } hkcal_f7args_t;

typedef struct hk_cal_info_t              /* digital to analog conversion info*/
         {
           int8_t param_code[MAXHKCALCHAR]; /* parameter code                    */           
           int8_t acronym[MAXHKCALCHAR];    /* acronym used by GSOC for parameter*/
           int8_t param_type[MAXHKCALCHAR]; /* parameter type:                   */ 
                                          /*         BAS = basic DER = derived */ 
           int32_t start_byte;               /* start byte within packet          */
                                          /* = 1 parameter does not depend on  */
                                          /*     any packet data               */
                                          /* != 1 parameter does depend on     */
                                          /*      packet data at this byte     */
           int32_t start_bit;                /* start bit for this parameter      */
           int32_t nbits;                    /* number of bits to be read         */
           int8_t int_type[MAXHKCALCHAR];   /* integer type:                     */
                                          /*     USG = unsigned ,SIG = signed */             
           int8_t output[MAXHKCALCHAR];     /* write output to file indicator    */
           int8_t output_format[MAXHKCALCHAR];/* output format I,F,E,T,H,B       */
           int32_t ndigits;                  /* number of digits after decimal    */
           int8_t cal_method[MAXHKCALCHAR]; /* calibration method                */
                                          /* scale(x)                          */
                                          /* poly(An,An-1,....,A0)             */
                                          /* intpol(x0,y0,x1,y1,...,xN,yN)     */
           hkcal_raw_t raw;               /* raw   function parameters         */
           hkcal_scale_t  scale;          /* scale function parameters         */
           hkcal_poly_t   poly;           /* poly  function parameters         */
           hkcal_f2args_t f2;             /* f2    function parameters         */
           hkcal_f3args_t f3;             /* f3    function parameters         */
           hkcal_f4args_t f4;             /* f4    function parameters         */
           hkcal_f5args_t f5;             /* f5    function parameters         */
           hkcal_f6args_t f6;             /* f6    function parameters         */
           hkcal_f7args_t f7;             /* f7    function parameters         */
         } hk_cal_info_t;

typedef struct thruster_info_t
         {
           double mass_dot;               /* nominal mass flow rate (kg/sec)   */
           double k_factor;               /* correction factor for mass_dot    */
           double mass_dot_sigma;         /* sigma for mass_dot     (kg/sec)   */
           double k_factor_sigma;         /* sigma for k_factor     (no dim)   */
           int32_t string_ndx;               /* index to which string thruster    */
                                          /* belongs 1 or 2. i.e sourc of gas  */
                                          /* tank 1 or 2                       */
           int32_t ntable;                   /* number of thruster firings listed */
                                          /* in input table. Only used for     */
                                          /* orbit thrusters 13+14             */
           double start_time[NTHRUST_MAX];/* start time thrust (GPS sec)       */
           double final_time[NTHRUST_MAX];/* stop  time thrust (GPS sec)       */
           double accum_time[NTHRUST_MAX];/* accumulate on time                */
           int32_t  thrust_count[NTHRUST_MAX];/* thrust event count               */
           double start_time_sigma[NTHRUST_MAX];/* sigma for start time (sec)  */
           double final_time_sigma[NTHRUST_MAX];/* sigma for stop  time (sec)  */
         } thruster_info_t;

typedef struct mas_com_info_t
         {
           int8_t sccsid[MAXHKCALCHAR];     /* sccs id of set up file            */
           int8_t grace_id;                 /* satellite id (A or B)             */
           double mass_zero;              /* initial satellite mass (no gas) kg*/
           double mass_gas_zero_1;        /* initial gas mass tank 1 in kg     */
           double pres_gas_zero_1;        /* initial pressure tank 1 in bar    */
           double temp_gas_zero_1;        /* initial temperature tank 1 (deg C)*/
           double mass_gas_zero_2;        /* initial gas mass tank 1 in kg     */
           double pres_gas_zero_2;        /* initial pressure tank 1 in bar    */
           double temp_gas_zero_2;        /* initial temperature tank 1 (deg C)*/
           double mass_zero_sigma;        /* sigma on mass_zero in kg          */
           double mass_gas_zero_1_sigma;  /* sigma on mass_gas_zero_1 in kg    */
           double pres_gas_zero_1_sigma;  /* sigma on pres_gas_zero_1 in bar   */
           double temp_gas_zero_1_sigma;  /* sigma on temp_gas_zero_1 in deg C */
           double mass_gas_zero_2_sigma;  /* sigma on mass_gas_zero_2 in kg    */
           double pres_gas_zero_2_sigma;  /* sigma on pres_gas_zero_2 in bar   */
           double temp_gas_zero_2_sigma;  /* sigma on temp_gas_zero_2 in deg C */
           thruster_info_t thruster[MAXTHRUSTERS]; /* info for all thrusters   */
         } mas_com_info_t;

typedef struct sensor_soe_info_t
         {
           int32_t   nevents;                /* number of events              */
           double *event_times;           /* pointer to event times array  */
           double *event_ndx;             /* pointer to event indices array*/
           int32_t   *nvalues;               /* number of values per event    */
           double **values;               /* values per event              */
         } sensor_soe_info_t;

TMpackets_t GRACE_TMpackets[MAXNRPACKETS];
mas_com_info_t GRACE_MAScom[2];

int32_t GRACE_TMpointer[MAXNRPACKETS];
int32_t GRACE_TMstats[MAXNRPACKETS][NRTMSTATS];

ICUDAcoef_t GRACE_ICUcal[NRICUSATS][NRICUBOARDS];

hk_cal_info_t GRACE_HKcal[MAXHKSAT][MAXHKCALLINES];

sensor_soe_info_t GRACE_SOE[2][NSOESENSORS];

/*>>>> prototypes <<<<*/

int32_t LoadTMproperties();
int32_t LoadHKCALproperties();
int32_t LoadSOEproperties(int8_t *soe_filename);
int32_t LoadMASCOMproperties();
int32_t LoadICUcal();
void PrintICUAP_info(FILE *dst, ICUAP_info_t *info,int32_t Informative,int32_t Flush);
int32_t EncodeICUpacket(uint8_t PktId,int8_t SuppPktId,uint16_t AP_Dlen,
                     int8_t *AP_Data,
                     int32_t RefTime, int8_t TS_AP_Suppid, ACC1A_t *screcs, 
                     int32_t *Nscrec, ACC1A_t *hkrecs, int32_t *Nhkrec,
                     int8_t GRACE_id,double Temperature,int32_t ICUbinary,
                     int32_t *NicuPkts, ICUAP_info_t *icupkt_info);
int32_t EncodeMAGpacket(uint8_t PktId,int8_t SuppPktId, uint16_t AP_Dlen,
                     int8_t *AP_Data, int32_t RefTime, int8_t TS_AP_Suppid, 
                     MAG1X_t *magrec, int8_t GRACE_id,double Temperature, 
                     int32_t MAGbinary);
int32_t EncodeAOCSpacket(uint8_t PktId,int8_t SuppPktId, uint16_t AP_Dlen,
                      int8_t *AP_Data, int32_t RefTime, int8_t TS_AP_Suppid,
                      THR1X_t *thrrec, int8_t GRACE_id,double Temperature,
                      int32_t THRbinary, double SCET_offset);
int32_t EncodeANARAWpacket(uint8_t PktId,int8_t SuppPktId, uint16_t AP_Dlen,
                        int8_t *AP_Data, int32_t RefTime, int8_t TS_AP_Suppid, 
                        TNK1X_t *tnkrecA, TNK1X_t *tnkrecB, int8_t GRACE_id,
                        double Temperature, int32_t TNKbinary, int32_t NewTNKrec);
int32_t EncodePWR_PCDUHKpacket(uint8_t PktId,int8_t SuppPktId,
                            uint16_t AP_Dlen,
                        int8_t *AP_Data, int32_t RefTime, int8_t TS_AP_Suppid,
                        TNK1X_t *tnkrecA, TNK1X_t *tnkrecB, int8_t GRACE_id,
                        double Temperature, int32_t TNKbinary, int32_t NewTNKrec);
int32_t EncodeHRTpacket(uint8_t PktId,int8_t SuppPktId, uint16_t AP_Dlen,
                     int8_t *AP_Data, int32_t RefTime, int8_t TS_AP_Suppid,
                     HRT1X_t *hrtrec, int8_t GRACE_id);
void FillTimStruct(int8_t GRACE_id, tm_sc_stat_t *tm_stat,TIM1X_t *tim);
int32_t PickLong(int8_t *buf);
int32_t Pick3Int(int8_t *buf);
int16_t Pick12bitInt(int8_t *buf);
int16_t PickShort(int8_t *buf);
int32_t ReadApplicationPacket(FILE *src,uint8_t *AP_PktId, int8_t *AP_SuppId,
                          uint16_t *AP_Dlen, int8_t *AP_Data);
int32_t WriteApplicationPacket(FILE *dst,uint8_t *AP_PktId, int8_t *AP_SuppId,
                          uint16_t *AP_Dlen, int8_t *AP_Data, int16_t byte_offset);
int32_t roundint(double x);
int32_t FindFirstGPSTag(FILE *src);
int32_t GetParamNdx(int8_t *name,int32_t isat);
int32_t InitHKComputationParameters(int8_t **Param_names, int32_t nparam_names,
                                 int32_t *var_ndx, int32_t *nvar_ndx, int32_t isat);
double raw(double raw_value);
double scale(double raw_value,double scale_factor);
double poly(double raw_value,double *coeff, int32_t ncoeff);
double f2(double raw_value,double arg1, double bias);
double f3(double raw_value,double arg1, double scale);
double f4(double raw_value,double arg1,double arg2, double arg3, double bias);
double f5(double arg1, double arg2);
double f6(double raw_value,double arg1,double arg2, double arg3);
double f7(double arg1, double arg2);
double ExtractValueFromHKTM(int32_t cal_ndx,int8_t *AP_data, int32_t isat);
int32_t GetThrAccumTimeCount(int32_t thr_number, int32_t ipsat, double T2000,
                          double *accum_time, int32_t *accum_count);
int32_t GetSOEParamNdx(int8_t *param_name);
void GetSOEParamName(int32_t isensor, int8_t *param_name);
int32_t GetActiveSensor(int8_t *sensor_name, int8_t *sat_name,
                     double start_time, double end_time,
                     int32_t *NactiveIntervals, double *start_time_interval,
                     double *end_time_interval,
                     double **soe_values, int32_t *nr_soe_values,
                     int32_t nint_max, int32_t nsoe_max);
double TankVolume(double pressure,int32_t TankId, int8_t GRACE_id);
double ComputeGN2GasMass(double presssure, double temperature, double Volume,
                         double He_percentage);
int32_t good_TS_packet(int8_t * pkt);
int32_t good_AP_packet(int8_t * pkt, int32_t Input_Dlen);
int32_t FindNextTimeGoodPkt(FILE *src,uint8_t *AP_PktId, int8_t *AP_SuppId,
                         uint16_t *AP_Dlen, int8_t *AP_Data,FILE *err,
                         uint16_t *ErrorPacketLength);


typedef struct {
                 int8_t parname[1000];
                 double *time;
                 double *rate;
                 int32_t n;
               } massflow_t;

typedef struct {
                 int32_t npar;
                 double mass_used;
                 double obs_weight;
                 double partials_read[MAXPARAM];
                 double reg_pressure;
                 double gas_temperature;
                 double thr_time;
                 double tank_time;
                 double mass_used_prev;
                 double tank_mass_calc;
                 int32_t nwords;
               } massflowreg_t;

int32_t ReadTankRegres(FILE *src,massflowreg_t *record);
int32_t WriteTankRegres(FILE *dst,massflowreg_t *record);
int32_t LoadMassFlowFromTdp(FILE *tdp, massflow_t *massflow, int8_t *tdp_filename);

#endif  /* _GRACEtelemetry_h_ */



