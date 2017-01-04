/* %Z% %M%      %I% %G%

   Purpose:	open file specified by environment variable and check for errors

   09/01/2000	Larry Romans	Created
*/
#include "GRACEdefs.h"
#include "GRACEprototype.h"

FILE *file_open(const int8_t *envvar, const int8_t *mode)
{
  static int8_t SccsId[] = "$Id: file_open.c,v 1.2 2009/06/06 22:28:26 glk Exp $";
  FILE *fp;
  int8_t *filename;
  filename = getenv(envvar);
  if (!filename) {
    fprintf(stderr, "Problem with environment variable %s (apparently not set)\n", envvar);
    exit(1);
  }
  fp = fopen(filename, mode);
  if (!fp) {
    fprintf(stderr,
      "Problem opening file %s with mode \"%s\" (from environment variable %s)\n",
      filename, mode, envvar);
    exit(1);
  }
  return(fp);
}

