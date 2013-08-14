
#include <stdio.h>
float* mat(const char * infile,int m,int n){
     int length = 0;
     float d;
     FILE    *fptr;
     /* Open the file */
     if ((fptr = fopen(infile,"r")) == (char)0) {
          fprintf(stderr,"Unable to open data file\n");
          printf("%s",infile);
          exit(0);
     }

     int row;
     int col;
     fscanf(fptr,"%i",&row);
     fscanf(fptr,"%i",&col);

     //printf("%i\n",row);
     //printf("%i\n\n\n",col);


     float* data = (float*) malloc(sizeof(float)*row*col);

     /* Read as many points as we can */
     while (fscanf(fptr,"%lf",&d) == 1) {
          if(length>row*col){
               printf("Malformed Data File\n");
          }
          data[length] = d;
          length++;
     }

     fclose(fptr);
     m = row;
     n = col;

}


int write(const char * outfile,int m,int n, float* data){

     int length = 0;
     float d;
     FILE    *fptr;
     /* Open the file */
     if ((fptr = fopen(outfile,"wt")) == NULL) {
          printf("Unable to open data file\n");
          printf("%s",outfile);
          return 1;
     }

        fprintf (fptr, "%i\n",m);
        fprintf (fptr, "%i\n",n);
        int i;
        for(i=0;i<m*n;i++)fprintf (fptr, "%f\n",data[i]);
        fclose (fptr);
        return 0;



}
