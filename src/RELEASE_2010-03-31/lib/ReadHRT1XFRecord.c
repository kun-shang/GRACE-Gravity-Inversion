#include "GRACEiolib.h"
#include "GRACEio_prototypes.h"

#define Success 1
#define Failure 0

static int8_t SccsId[] = "$Id: ReadHRT1XFRecord.c,v 1.3 2009/06/06 22:15:38 glk Exp $";

boolean ReadHRT1XFRecord(FILE *src, HRT1X_t *record)
/*----------------------------------------------------------------------------->
/ purpose: read the High Resolution Temperature data from flee 
/          pointer src
/
/ Initial version: Gerhard L.H. Kruizinga      08/05/04
/
/
/ input:  *src    Pointer to HRT1X Data Format File
/ output: *record Pointer to HRT1X Data struct (HRT1X_t)
/
/ return:      1       normal return
/              0       End Of File reached
<-----------------------------------------------------------------------------*/
{
  uint32_t time_intg;

  int32_t    retrn,i;

/*----------------------------------------------------------------------------->
/ Test for EOF
<-----------------------------------------------------------------------------*/

  retrn = fread_grace(&time_intg,sizeof(record->time_intg),1,src);
  if (feof(src) != 0) return False;
  if (retrn != 1)
  {
    printf("HRT1X: Error reading field 'time_intg' \n");
    return Failure;
  }
  record->time_intg = time_intg;
  
/*----------------------------------------------------------------------------->
/ Read remaining record elements from src
<-----------------------------------------------------------------------------*/

  if (fread_grace(&record->time_frac,sizeof(record->time_frac),1,src) != 1)
  {
    printf("HRT1X: Error reading field 'time_frac' \n");
    return Failure;
  }

  if (fread_grace(&record->time_ref,sizeof(record->time_ref),1,src) != 1)
  {
    printf("HRT1X: Error reading field 'time_ref' \n");
    return Failure;
  }

  if (fread_grace(&record->GRACE_id,sizeof(record->GRACE_id),1,src) != 1)
  {
    printf("HRT1X: Error reading field 'GRACE_id' \n");
    return Failure;
  }

  if (fread_grace(&record->TEMP_MEP_neg_y,sizeof(record->TEMP_MEP_neg_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MEP_neg_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MEP_pos_y,sizeof(record->TEMP_MEP_pos_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MEP_pos_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MEPm,sizeof(record->TEMP_MEPm),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MEPm \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ICU,sizeof(record->TEMP_ICU),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ICU \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ICU_red,sizeof(record->TEMP_ICU_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ICU_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ACC_neg_z,sizeof(record->TEMP_ACC_neg_z),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ACC_neg_z \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ACC_pos_z,sizeof(record->TEMP_ACC_pos_z),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ACC_pos_z \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_pos_x,sizeof(record->TEMP_CFRP_pos_x),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_pos_x \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_pos_x_red,sizeof(record->TEMP_CFRP_pos_x_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_pos_x_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_neg_x,sizeof(record->TEMP_CFRP_neg_x),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_neg_x \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_neg_x_red,sizeof(record->TEMP_CFRP_neg_x_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_neg_x_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_neg_y,sizeof(record->TEMP_CFRP_neg_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_neg_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_CFRP_neg_y_red,sizeof(record->TEMP_CFRP_neg_y_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_CFRP_neg_y_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ACCSen,sizeof(record->TEMP_ACCSen),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ACCSen \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_ICU_spec,sizeof(record->TEMP_ICU_spec),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_ICU_spec \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MWA_neg_y,sizeof(record->TEMP_MWA_neg_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MWA_neg_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MWA_neg_yoff,sizeof(record->TEMP_MWA_neg_yoff),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MWA_neg_yoff \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MWA_pos_y,sizeof(record->TEMP_MWA_pos_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MWA_pos_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_MWA_pos_yoff,sizeof(record->TEMP_MWA_pos_yoff),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_MWA_pos_yoff \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_Horn_pos_x,sizeof(record->TEMP_Horn_pos_x),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_Horn_pos_x \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_Horn_pos_x_red,sizeof(record->TEMP_Horn_pos_x_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_Horn_pos_x_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_HornPl,sizeof(record->TEMP_HornPl),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_HornPl \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_HornPl_red,sizeof(record->TEMP_HornPl_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_HornPl_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_HMWA_neg_y,sizeof(record->TEMP_HMWA_neg_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_HMWA_neg_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_HMWA_pos_y,sizeof(record->TEMP_HMWA_pos_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_HMWA_pos_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_RFSamp,sizeof(record->TEMP_RFSamp),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_RFSamp \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_USO_neg_y,sizeof(record->TEMP_USO_neg_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_USO_neg_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_USO_neg_y_red,sizeof(record->TEMP_USO_neg_y_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_USO_neg_y_red \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_USO_pos_y,sizeof(record->TEMP_USO_pos_y),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_USO_pos_y \n" ) ; 
    return Failure ; 
  }

  if (fread_grace(&record->TEMP_USO_pos_y_red,sizeof(record->TEMP_USO_pos_y_red),1,src) != 1)
  {
    printf("HRT1X: Error reading field TEMP_USO_pos_y_red \n" ) ; 
    return Failure ; 
  }


  if (fread_grace(&record->qualflg,sizeof(record->qualflg),1,src) != 1)
  {
     printf("HRT1X: Error reading field 'qualflg' \n");
     return Failure;
  }

  return True;

}
