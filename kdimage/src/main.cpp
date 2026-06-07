#include "kdimage.h"


int main(int argc, char *argv[])
{

    if(argc!=2)
    {
        fprintf(stderr,"kdimage <data.dat>\n");
        return 0;
    }

    KDIMAGE kd;

    fprintf(stderr,"Data file: %s\n",argv[1]);

    FILE *f=fopen(argv[1],"rt");

    kdimage_init(f,&kd);
    kdimage_print(&kd);
    kdimage_free(&kd);

    fclose(f);

    return 0;

}
