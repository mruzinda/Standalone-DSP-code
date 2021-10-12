#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#define NUM_COMMANDS 2

// This script is used to test different plotting methods for use in coherent beamformer output analysis

int main(){
  ///*
  int test_flag = 1; // What kind of plot do you want?
  if (test_flag == 0){ // Circle plot
      char * commandsForGnuplot[] = {"set title \"Plot of a circle\"", "plot 'circle.txt'"};
	  FILE *fp=NULL;
	  fp=fopen("circle.txt","w");
	  double r;
      double x,y,x0,y0;
      printf("Enter the radius of the circle to be plotted: ");
      scanf("%lf",&r);
      printf("Enter the x and y-coordinates of the center: ");
      scanf("%lf%lf",&x0,&y0);
	  
      for(y=y0-r;y<=y0+r;y=y+0.1){
        x=sqrt(r*r-(y-y0)*(y-y0))+x0; 
        fprintf(fp,"%lf\t %lf\n",x,y);
      }
      for(y=y0+r;y>=y0-r;y=y-0.1){
        x=-sqrt(r*r-(y-y0)*(y-y0))+x0; 
        fprintf(fp,"%lf\t %lf\n",x,y);
      }
	  FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
	  
	  for(int i=0; i<NUM_COMMANDS; i++){
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
	  }
  
      fflush(gnuplotPipe);
  }
  if(test_flag == 1){ // Intensity or heat map
      char * commandsForGnuplot[] = {"set title \"A heat map\"", "plot 'heatmap.txt' matrix with image"};
	  FILE *fp=NULL;
	  fp=fopen("heatmap.txt","w");
      double x,y,z;
	  
      for(z=0;z<3;z++){
	    for(y=0;y<3;y++){
		  for(x=0;x<3;x++){
		    fprintf(fp,"%lf\t %lf\t %lf\n",x,y,z);
	      }
        }
	  }
  
      FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
  
      //fprintf(gnuplotPipe, "plot 'heatmap.txt' matrix with image \n");
	  
	  for(int i=0; i<NUM_COMMANDS; i++){
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
	  }
  
      fflush(gnuplotPipe);
  }
  if(test_flag == 2){ // 3D plot of a sphere
	  char * commandsForGnuplot[] = {"set title \"Plot of a sphere\"", "splot 'circle_map.txt' w 1"};
	  FILE *fp=NULL;
      fp=fopen("circle_map.txt","w");
      double r;
      double x,y,z,x0,y0,z0;
      printf("Enter the radius of the circle to be plotted: ");
      scanf("%lf",&r);
      printf("Enter the x, y, and z-coordinates of the center: ");
      scanf("%lf%lf%lf",&x0,&y0,&z0);
      for(y=y0-r;y<=y0+r;y=y+0.1){
	    for(z=z0-r;z<=z0+r;z=z+0.1){
		  x=sqrt(r*r-(y-y0)*(y-y0)-(z-z0)*(z-z0))+x0; 
		  fprintf(fp,"%lf\t %lf\t %lf\n",x,y,z);
	    }
      }
      for(y=y0+r;y>=y0-r;y=y-0.1){
	    for(z=z0+r;z>=z0-r;z=z-0.1){
		  x=-sqrt(r*r-(y-y0)*(y-y0)-(z-z0)*(z-z0))+x0; 
		  fprintf(fp,"%lf\t %lf\t %lf\n",x,y,z);
	    }
      }
  
      FILE * gnuplotPipe = popen ("gnuplot -persistent", "w");
  
      //fprintf(gnuplotPipe, "splot 'circle_map.txt' w l \n");
	  
	  for(int i=0; i<NUM_COMMANDS; i++){
        fprintf(gnuplotPipe, "%s \n", commandsForGnuplot[i]);
	  }
  
      fflush(gnuplotPipe);
  
  }
  
  return 0;
}