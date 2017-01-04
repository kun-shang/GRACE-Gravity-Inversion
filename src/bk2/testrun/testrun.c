#include <stdio.h>
#include <math.h>
#include "coord.h"
#include "novas.h"

int main()
{

    double tjd[2], JD0, lib[6];
    int i;

    long int target, center;


    printf ("test begin\n");


    InfStruct info;


    target = 15;
    center = 0;

//    dpleph_(jd, &target, &center, lib);



    JD0 = julian_date (2013, 5, 1,0);

    eop_open("/0/home/shang.34/inputs/eopc04_08_IAU2000.62-now",JD0- 2400000.5, JD0 + 30 - 2400000.5);


    for (i = 0; i <= 86400; i=i+3600)
    {

        tjd[0] = JD0;    tjd[1] = (i + 67.184) / 86400.0;

        getinfo(tjd, 3, &info);

        printf ("%f\t%.20f\n", info.utc, info.c_ie[0]);


    }

    for (i = 0; i <= 86400*30; i=i+86400)
    {

        tjd[0] = JD0;    tjd[1] = (i + 67.184) / 86400.0;
        dpleph_(tjd, &target, &center, lib);


        getinfo(tjd, 9, &info);

        printf ("%f\t%.20f\n", info.utc/86400, info.c_ie[0]);

    }

    eop_close();


    return 0;
}
