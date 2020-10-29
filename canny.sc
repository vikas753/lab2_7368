/*******************************************************************************
* --------------------------------------------
*(c) 2001 University of South Florida, Tampa
* Use, or copying without permission prohibited.
* PERMISSION TO USE
* In transmitting this software, permission to use for research and
* educational purposes is hereby granted.  This software may be copied for
* archival and backup purposes only.  This software may not be transmitted
* to a third party without prior permission of the copyright holder. This
* permission may be granted only by Mike Heath or Prof. Sudeep Sarkar of
* University of South Florida (sarkar@csee.usf.edu). Acknowledgment as
* appropriate is respectfully requested.
* 
*  Heath, M., Sarkar, S., Sanocki, T., and Bowyer, K. Comparison of edge
*    detectors: a methodology and initial study, Computer Vision and Image
*    Understanding 69 (1), 38-54, January 1998.
*  Heath, M., Sarkar, S., Sanocki, T. and Bowyer, K.W. A Robust Visual
*    Method for Assessing the Relative Performance of Edge Detection
*    Algorithms, IEEE Transactions on Pattern Analysis and Machine
*    Intelligence 19 (12),  1338-1359, December 1997.
*  ------------------------------------------------------
*
* PROGRAM: canny_edge
* PURPOSE: This program implements a "Canny" edge detector. The processing
* steps are as follows:
*
*   1) Convolve the image with a separable gaussian filter.
*   2) Take the dx and dy the first derivatives using [-1,0,1] and [1,0,-1]'.
*   3) Compute the magnitude: sqrt(dx*dx+dy*dy).
*   4) Perform non-maximal suppression.
*   5) Perform hysteresis.
*
* The user must input three parameters. These are as follows:
*
*   sigma = The standard deviation of the gaussian smoothing filter.
*   tlow  = Specifies the low value to use in hysteresis. This is a 
*           fraction (0-1) of the computed high threshold edge strength value.
*   thigh = Specifies the high value to use in hysteresis. This fraction (0-1)
*           specifies the percentage point in a histogram of the gradient of
*           the magnitude. Magnitude values of zero are not counted in the
*           histogram.
*
* NAME: Mike Heath
*       Computer Vision Laboratory
*       University of South Floeida
*       heath@csee.usf.edu
*
* DATE: 2/15/96
*
* Modified: 5/17/96 - To write out a floating point RAW headerless file of
*                     the edge gradient "up the edge" where the angle is
*                     defined in radians counterclockwise from the x direction.
*                     (Mike Heath)
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <c_typed_queue.sh>

#define Q_SIZE 2ul


#define VERBOSE 0
#define BOOSTBLURFACTOR 90.0


#define MAX_ROWS 240
#define MAX_COLS 320

#define SIGMA 0.6f
#define TLOW  0.3f
#define THIGH 0.8f

// Data structure for holding image contents
typedef struct image_data
{
  int rows;
  int cols;
  int terminate_value;
  unsigned char img[MAX_ROWS*MAX_COLS];
}image_data_t;


// Data structure for holding smoothened or derivate of image contents
typedef struct image_data_short
{
  short img[MAX_ROWS*MAX_COLS];
}image_data_short_t;


// define interface and implementation of typed queue
DEFINE_IC_TYPED_QUEUE(imagedata, image_data_t)


#define WINDOW_SIZE 21

unsigned int read_pgm_image(char *infilename, unsigned char **image, int *rows,
    int *cols,image_data_t* read_image_buffer);
int write_pgm_image(char *outfilename, unsigned char *image, int rows,
    int cols, char *comment, int maxval);

void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
        float tlow, float thigh, unsigned char *edge);
void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,unsigned char *result); 

#define TERMINATE_SENTINEL 0x5a5a5a5a

#define MAX_ARGS 10

// These are the only set of global variables which will be used to set 
// system arguments of stimulus instance
int argc_glo;
char *argv_glo[MAX_ARGS];
  
// sender behavior with input image
behavior Stimulus(i_imagedata_sender ch) {
  image_data_t stimulus_beh_memory;
  int argc,i=0;
  char *argv[MAX_ARGS];
  
  void main(void)
  {
     char *infilename = 0;  /* Name of the input image */
     unsigned char *image;     /* The input image */
     int iter_index = 1;           /* The dimensions of the image. */

     
     argc = argc_glo;
     for(;i<argc;i++)
     {
       argv[i] = argv_glo[i];
     }
     
     /****************************************************************************
     * Get the command line arguments.
     ****************************************************************************/
     if(argc < 2){
       fprintf(stderr,"\n<USAGE> %s image \n",argv[0]);
       fprintf(stderr,"\n      image:      An image to process. Must be in ");
       fprintf(stderr,"PGM format.\n");
       exit(1);
     }
     
     while(iter_index<argc)
     {
       infilename = argv[iter_index];
   
       /****************************************************************************
       * Read in the image. This read function allocates memory for the image.
       ****************************************************************************/
       if(VERBOSE) printf("Reading the image %s.\n", infilename);
       if(read_pgm_image(infilename, &image, &stimulus_beh_memory.rows, &stimulus_beh_memory.cols,&stimulus_beh_memory) == 0){
         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
         exit(1);
       }
       
       stimulus_beh_memory.terminate_value = 0;
       ch.send(stimulus_beh_memory);
       iter_index++;
     }
     stimulus_beh_memory.terminate_value = TERMINATE_SENTINEL;
     ch.send(stimulus_beh_memory); 
   }
};

// Data input block which takes in an image as an input
behavior Datain(i_imagedata_receiver rec_ch,i_imagedata_sender send_ch)
{
  image_data_t datain_beh_memory;
  void main(void)
  {
    while(1)
    {
      rec_ch.receive(&datain_beh_memory);
      send_ch.send(datain_beh_memory);  
      if(datain_beh_memory.terminate_value == TERMINATE_SENTINEL)
      {
        break;
      }
      
    }
  }  
};


// Gaussian Smooth : child behavior , which smoothens an image and produces a smoothened output
behavior gaussian_smooth(in image_data_t image_arg,out image_data_short_t smoothedim)
{
  int i;
  float x, fx;
  
  int r, c, rr, cc,     /* Counter variables. */
      windowsize,        /* Dimension of the gaussian kernel. */
      center;            /* Half of the windowsize. */
      
  float *tempim,        /* Buffer for separable filter gaussian smoothing. */
         *kernel,        /* A one dimensional gaussian kernel. */
         dot,            /* Dot product summing variable. */
         sum,sigma;            /* Sum of the kernel weights variable. */
         
  int rows,cols;
  unsigned char* image;
  float tempim_static_buffer[MAX_ROWS*MAX_COLS];
  float kernel_static_buffer[WINDOW_SIZE];
   
  void main(void)
  {
   image = image_arg.img;
   
   tempim = tempim_static_buffer;
   kernel = kernel_static_buffer;
   
   rows = MAX_ROWS;
   cols = MAX_COLS;
   sigma = 0.6;
   

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the gaussian smoothing kernel.\n");

   windowsize = 1 + 2 * ceil(2.5 * sigma);
   center = windowsize / 2;

   if(VERBOSE) printf("      The kernel has %d elements.\n", windowsize);
   
   for(i=0;i<windowsize;i++){
      x = (float)(i - center);
      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
      kernel[i] = fx;
      sum += fx;
   }

   for(i=0;i<windowsize;i++) kernel[i] /= sum;

   if(VERBOSE){
      printf("The filter coefficients are:\n");
      for(i=0;i<windowsize;i++)
         printf("kernel[%d] = %f\n", i, kernel[i]);
   }   center = windowsize / 2;

   /****************************************************************************
   * Assign buffer image and the smoothed image static buffers.
   ****************************************************************************/
   
   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
   for(r=0;r<rows;r++){
      for(c=0;c<cols;c++){
         dot = 0.0;
         sum = 0.0;
         for(cc=(-center);cc<=center;cc++){
            if(((c+cc) >= 0) && ((c+cc) < cols)){
               dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
               sum += kernel[center+cc];
            }
         }
         tempim[r*cols+c] = dot/sum;
      }
   }

   /****************************************************************************
   * Blur in the y - direction.
   ****************************************************************************/
   if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
   for(c=0;c<cols;c++){
      for(r=0;r<rows;r++){
         sum = 0.0;
         dot = 0.0;
         for(rr=(-center);rr<=center;rr++){
            if(((r+rr) >= 0) && ((r+rr) < rows)){
               dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
               sum += kernel[center+rr];
            }
         }
         smoothedim.img[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
      }
   }
    
  }    
}; 

// Derivative x,y behavior : Performs a derivative on x and y-axis of an image and sends out on a port.
behavior Derivative_x_y(in image_data_short_t smooth_data_arg,out image_data_short_t delta_x,out image_data_short_t delta_y)
{
  int r, c, pos, rows, cols;
  short int *smoothedim;
  
  void main(void)
  {
   smoothedim = smooth_data_arg.img;
   
   /****************************************************************************
   * Assign static buffers to store the derivatives.
   ****************************************************************************/
   rows = MAX_ROWS;
   cols = MAX_COLS;
   
   /****************************************************************************
   * Compute the x-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the X-direction derivative.\n");
   for(r=0;r<rows;r++){
      pos = r * cols;
      delta_x.img[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      for(c=1;c<(cols-1);c++,pos++){
         delta_x.img[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      }
      delta_x.img[pos] = smoothedim[pos] - smoothedim[pos-1];
   }

   /****************************************************************************
   * Compute the y-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
   for(c=0;c<cols;c++){
      pos = c;
      delta_y.img[pos] = smoothedim[pos+cols] - smoothedim[pos];
      pos += cols;
      for(r=1;r<(rows-1);r++,pos+=cols){
         delta_y.img[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
      }
      delta_y.img[pos] = smoothedim[pos] - smoothedim[pos-cols];
   }    
  }
};

// Magnitude x,y : Computes the magnitude in X and Y direction and sends output on a port
behavior Magnitude_x_y(in image_data_short_t delta_x_img, in image_data_short_t delta_y_img,out image_data_short_t magnitude_img)
{
  int r, c, pos, sq1, sq2,rows ,cols;
  short int *delta_x; 
  short int *delta_y;
  
  void main(void)
  {
    delta_x = delta_x_img.img;
    delta_y = delta_y_img.img;
    /****************************************************************************
    * Assign an image to store the magnitude of the gradient.
    ****************************************************************************/
    rows = MAX_ROWS;
    cols = MAX_COLS;

    for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
         magnitude_img.img[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
      }
    }

  }  
};

// Non maximum suppression : Performs a non maximal suppression on an image alongwith x and y gradients and sends output on a port
behavior Non_max_supp(in image_data_short_t mag_img,in image_data_short_t gradx_img,in image_data_short_t grady_img,out image_data_t result_img)
{
  int rows,cols;
  unsigned char* result_arg;
  unsigned char nms_static_buffer_temp[MAX_ROWS*MAX_COLS];
  int i=0;
  short *mag;
  short *gradx;
  short *grady;
  
  void main(void)
  {
    rows = MAX_ROWS;
    cols = MAX_COLS;
    mag = mag_img.img;
    gradx = gradx_img.img;
    grady = grady_img.img;
    
    non_max_supp(mag, gradx, grady, rows, cols, nms_static_buffer_temp);
    for(i=0;i<(MAX_ROWS*MAX_COLS);i++)
    {
      result_img.img[i] = nms_static_buffer_temp[i];    
    } 
  }    
}; 

// Apply Hysterisi behavior : Computes hysterisis on an image and sends output on a port
behavior Apply_hysteresis(in image_data_short_t mag_img, in image_data_t nms_img,out image_data_t edge_img)
{
  int rows,cols,i;
  float tlow,thigh;
  
  short int *mag; 
  unsigned char *nms;
  
  unsigned char hysterisis_static_buffer_temp[MAX_ROWS*MAX_COLS];
  
  void main(void)
  {
    mag = mag_img.img;
    nms = nms_img.img;
    
    rows = MAX_ROWS;
    cols = MAX_COLS;
    tlow = TLOW;
    thigh = THIGH;
    
    apply_hysteresis(mag, nms, rows, cols, tlow, thigh, hysterisis_static_buffer_temp);
    for(i=0;i<(MAX_ROWS*MAX_COLS);i++)
    {
      edge_img.img[i] = hysterisis_static_buffer_temp[i];    
    } 
  }    
};

// Device Under test which runs canny algorithm
behavior DUT(i_imagedata_receiver rec_ch,i_imagedata_sender send_ch)
{
  image_data_t dut_beh_memory_in,dut_beh_memory_out;

  image_data_short_t smoothedim_data;
  
  image_data_short_t deltaX_data;
  image_data_short_t deltaY_data;  
  image_data_short_t magnitude_data;  
  image_data_t       nms_data;
  
  gaussian_smooth   gaussian_smooth_instance(dut_beh_memory_in,smoothedim_data);
  Derivative_x_y    derivative_x_y_instance(smoothedim_data,deltaX_data,deltaY_data);
  Magnitude_x_y     magnitude_x_y_instance(deltaX_data,deltaY_data,magnitude_data);
  Non_max_supp      non_max_supp_instance(magnitude_data,deltaX_data,deltaY_data,nms_data);
  Apply_hysteresis  apply_hysterisis_instance(magnitude_data,nms_data,dut_beh_memory_out);
                                   
  void main(void)
  {

    while(1)
    {
      rec_ch.receive(&dut_beh_memory_in);
      
      dut_beh_memory_out.rows = dut_beh_memory_in.rows;
      dut_beh_memory_out.cols = dut_beh_memory_in.cols;
      
      if(dut_beh_memory_in.terminate_value == TERMINATE_SENTINEL)
      {
        dut_beh_memory_out.terminate_value = TERMINATE_SENTINEL;
        send_ch.send(dut_beh_memory_out);
        break;
      }
      
      // Sequential execution of all child behaviors of DUT
      gaussian_smooth_instance;
      derivative_x_y_instance;
      magnitude_x_y_instance;
      non_max_supp_instance;
      apply_hysterisis_instance;
      
      send_ch.send(dut_beh_memory_out);
        
    }
  }  
};

// Dataout block of platform to send image to monitor
behavior Dataout(i_imagedata_receiver rec_ch,i_imagedata_sender send_ch)
{
  image_data_t dataout_beh_memory;
  void main(void)
  {
    while(1)
    {
      rec_ch.receive(&dataout_beh_memory);
      
      send_ch.send(dataout_beh_memory);  
      if(dataout_beh_memory.terminate_value == TERMINATE_SENTINEL)
      {
        break;
      }
    
      
    }
  }  
};

// Monitor behavior , for writing the output
behavior Monitor(i_imagedata_receiver rec_ch)
{
  int image_index = 0;
  image_data_t monitor_beh_memory;
  char blank_comment[128]; /* Constant comment */
  char outfilename[128];    /* Name of the output "edge" image */

  void main(void)
  {
    while(1)
    {
      rec_ch.receive(&monitor_beh_memory);
      
      if(monitor_beh_memory.terminate_value == TERMINATE_SENTINEL)
      {
        break;
      }
      
      sprintf(outfilename, "sldl_output_%d_sc.pgm", image_index);
      blank_comment[0] = 0;
      image_index++;
      if(VERBOSE) printf("Writing the edge iname in the file %s. rows : %d , cols : %d \n", outfilename , monitor_beh_memory.rows, monitor_beh_memory.cols);
        if(write_pgm_image(outfilename, &monitor_beh_memory.img[0], monitor_beh_memory.rows, monitor_beh_memory.cols, blank_comment , 255) == 0){
        fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
        exit(1);
      }        
    }
  }  
};

// Platform behavior , spawns datain , dut and dataout
behavior Platform(i_imagedata_receiver rec_i,i_imagedata_sender send_i)
{
	c_imagedata_queue datain_dut_queue(Q_SIZE),dut_dataout_queue(Q_SIZE);
 
	Datain datain_instance(rec_i,datain_dut_queue);
  DUT dut_instance(datain_dut_queue,dut_dataout_queue);
  Dataout dataout_instance(dut_dataout_queue,send_i);
  
	void main(void)	{
    par
    {
			datain_instance;
			dut_instance;
      dataout_instance;
    }
	}
};



// main behavior ( spawns stimulus , platform and monitor )
behavior Main()
{
	c_imagedata_queue stimulus_platform_queue(Q_SIZE),platform_monitor_queue(Q_SIZE);
	int i=0;
	Platform platform_instance(stimulus_platform_queue,platform_monitor_queue);
  Monitor  monitor_instance(platform_monitor_queue);
	Stimulus stimulus_instance(stimulus_platform_queue);
 
  int main(int argc, char *argv[])	{
    // Workaround for setting system arguments onto global variables and then using them in stimulus
    argc_glo = argc;
    for(;i<argc;i++)
    {
      argv_glo[i] = argv[i];
    }
    
		par {
		  stimulus_instance;
	    platform_instance;
      monitor_instance;
		}
		return 0;
	}
};

/*******************************************************************************
* FILE: hysteresis.c
* This code was re-written by Mike Heath from original code obtained indirectly
* from Michigan State University. heath@csee.usf.edu (Re-written in 1996).
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0

/*******************************************************************************
* PROCEDURE: follow_edges
* PURPOSE: This procedure edges is a recursive routine that traces edgs along
* all paths whose magnitude values remain above some specifyable lower
* threshhold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
   int cols)
{
   short *tempmagptr;
   unsigned char *tempmapptr;
   int i;
   int x[8] = {1,1,0,-1,-1,-1,0,1},
       y[8] = {0,1,1,1,0,-1,-1,-1};

   for(i=0;i<8;i++){
      tempmapptr = edgemapptr - y[i]*cols + x[i];
      tempmagptr = edgemagptr - y[i]*cols + x[i];

      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
         *tempmapptr = (unsigned char) EDGE;
         follow_edges(tempmapptr,tempmagptr, lowval, cols);
      }
   }
}

/*******************************************************************************
* PROCEDURE: apply_hysteresis
* PURPOSE: This routine finds edges that are above some high threshhold or
* are connected to a high pixel by a path of pixels greater than a low
* threshold.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
	float tlow, float thigh, unsigned char *edge)
{
   int r, c, pos, numedges, highcount, lowthreshold, highthreshold,hist[32768];
   short int maximum_mag;

   /****************************************************************************
   * Initialize the edge map to possible edges everywhere the non-maximal
   * suppression suggested there could be an edge except for the border. At
   * the border we say there can not be an edge because it makes the
   * follow_edges algorithm more efficient to not worry about tracking an
   * edge off the side of the image.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
	 else edge[pos] = NOEDGE;
      }
   }

   for(r=0,pos=0;r<rows;r++,pos+=cols){
      edge[pos] = NOEDGE;
      edge[pos+cols-1] = NOEDGE;
   }
   pos = (rows-1) * cols;
   for(c=0;c<cols;c++,pos++){
      edge[c] = NOEDGE;
      edge[pos] = NOEDGE;
   }

   /****************************************************************************
   * Compute the histogram of the magnitude image. Then use the histogram to
   * compute hysteresis thresholds.
   ****************************************************************************/
   for(r=0;r<32768;r++) hist[r] = 0;
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
      }
   }

   /****************************************************************************
   * Compute the number of pixels that passed the nonmaximal suppression.
   ****************************************************************************/
   for(r=1,numedges=0;r<32768;r++){
      if(hist[r] != 0) maximum_mag = r;
      numedges += hist[r];
   }

   highcount = (int)(numedges * thigh + 0.5);

   /****************************************************************************
   * Compute the high threshold value as the (100 * thigh) percentage point
   * in the magnitude of the gradient histogram of all the pixels that passes
   * non-maximal suppression. Then calculate the low threshold as a fraction
   * of the computed high threshold value. John Canny said in his paper
   * "A Computational Approach to Edge Detection" that "The ratio of the
   * high to low threshold in the implementation is in the range two or three
   * to one." That means that in terms of this implementation, we should
   * choose tlow ~= 0.5 or 0.33333.
   ****************************************************************************/
   r = 1;
   numedges = hist[1];
   while((r<(maximum_mag-1)) && (numedges < highcount)){
      r++;
      numedges += hist[r];
   }
   highthreshold = r;
   lowthreshold = (int)(highthreshold * tlow + 0.5);

   if(VERBOSE){
      printf("The input low and high fractions of %f and %f computed to\n",
	 tlow, thigh);
      printf("magnitude of the gradient threshold values of: %d %d\n",
	 lowthreshold, highthreshold);
   }

   /****************************************************************************
   * This loop looks for pixels above the highthreshold to locate edges and
   * then calls follow_edges to continue the edge.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++){
	 if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
            edge[pos] = EDGE;
            follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
	 }
      }
   }

   /****************************************************************************
   * Set all the remaining possible edges to non-edges.
   ****************************************************************************/
   for(r=0,pos=0;r<rows;r++){
      for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
   }
}

/*******************************************************************************
* PROCEDURE: non_max_supp
* PURPOSE: This routine applies non-maximal suppression to the magnitude of
* the gradient image.
* NAME: Mike Heath
* DATE: 2/15/96
*******************************************************************************/
void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,unsigned char *result) 
{
    int rowcount, colcount,count;
    short *magrowptr,*magptr;
    short *gxrowptr,*gxptr;
    short *gyrowptr,*gyptr,z1,z2;
    short m00,gx,gy;
    float mag1,mag2,xperp,yperp;
    unsigned char *resultrowptr, *resultptr;
    

   /****************************************************************************
   * Zero the edges of the result image.
   ****************************************************************************/
    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1); 
        count<ncols; resultptr++,resultrowptr++,count++){
        *resultrowptr = *resultptr = (unsigned char) 0;
    }

    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
        *resultptr = *resultrowptr = (unsigned char) 0;
    }

   /****************************************************************************
   * Suppress non-maximum points.
   ****************************************************************************/
   for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
      gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
      rowcount<=nrows-2; 
      rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
      resultrowptr+=ncols){   
      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
         resultptr=resultrowptr;colcount<=ncols-2; 
         colcount++,magptr++,gxptr++,gyptr++,resultptr++){   
         m00 = *magptr;
         if(m00 == 0){
            *resultptr = (unsigned char) NOEDGE;
         }
         else{
            xperp = -(gx = *gxptr)/((float)m00);
            yperp = (gy = *gyptr)/((float)m00);
         }

         if(gx >= 0){
            if(gy >= 0){
                    if (gx >= gy)
                    {  
                        /* 111 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                        
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {    
                        /* 110 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp; 
                    }
                }
                else
                {
                    if (gx >= -gy)
                    {
                        /* 101 */
                        /* Left point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;
            
                        /* Right point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {    
                        /* 100 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp; 
                    }
                }
            }
            else
            {
                if ((gy = *gyptr) >= 0)
                {
                    if (-gx >= gy)
                    {          
                        /* 011 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
                    }
                    else
                    {
                        /* 010 */
                        /* Left point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols + 1);

                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                        /* Right point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols - 1);

                        mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
                    }
                }
                else
                {
                    if (-gx > -gy)
                    {
                        /* 001 */
                        /* Left point */
                        z1 = *(magptr + 1);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                        /* Right point */
                        z1 = *(magptr - 1);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
                    }
                    else
                    {
                        /* 000 */
                        /* Left point */
                        z1 = *(magptr + ncols);
                        z2 = *(magptr + ncols + 1);

                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                        /* Right point */
                        z1 = *(magptr - ncols);
                        z2 = *(magptr - ncols - 1);

                        mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
                    }
                }
            } 

            /* Now determine if the current point is a maximum point */

            if ((mag1 > 0.0) || (mag2 > 0.0))
            {
                *resultptr = (unsigned char) NOEDGE;
            }
            else
            {    
                if (mag2 == 0.0)
                    *resultptr = (unsigned char) NOEDGE;
                else
                    *resultptr = (unsigned char) POSSIBLE_EDGE;
            }
        } 
    }
}
/*******************************************************************************
* FILE: pgm_io.c
* This code was written by Mike Heath. heath@csee.usf.edu (in 1995).
*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/******************************************************************************
* Function: read_pgm_image
* Purpose: This function reads in an image in PGM format. The image can be
* read in from either a file or from standard input. The image is only read
* from standard input when infilename = NULL. Because the PGM format includes
* the number of columns and the number of rows in the image, these are read
* from the file. Memory to store the image is allocated in this function.
* All comments in the header are discarded in the process of reading the
* image. Upon failure, this function returns 0, upon sucess it returns 1.
******************************************************************************/
unsigned int read_pgm_image(char *infilename, unsigned char **image, int *rows,
    int *cols,image_data_t* read_image_buffer)
{
   FILE *fp;
   char buf[71];

   /***************************************************************************
   * Open the input image file for reading if a filename was given. If no
   * filename was provided, set fp to read from standard input.
   ***************************************************************************/
   if(infilename == 0) fp = stdin;
   else{
      if((fp = fopen(infilename, "r")) == 0){
         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
            infilename);
         return(0);
      }
   }

   /***************************************************************************
   * Verify that the image is in PGM format, read in the number of columns
   * and rows in the image and scan past all of the header information.
   ***************************************************************************/
   fgets(buf, 70, fp);
   if(strncmp(buf,"P5",2) != 0){
      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
      fprintf(stderr, "read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   sscanf(buf, "%d %d", cols, rows);
   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
   
   if(((*rows) > MAX_ROWS) || ((*cols) > MAX_COLS))
   {
     fprintf(stderr, "Image size exceeded in read_pgm_image().\n ");
   }
   
   /***************************************************************************
   * Assign the max buffer to store image
   ***************************************************************************/

   (*image) = (unsigned char *) read_image_buffer->img;

   if((unsigned int)(*rows) != fread((*image), (*cols), (*rows), fp)){
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      if(fp != stdin) fclose(fp);
      return(0);
   }

   if(fp != stdin) fclose(fp);
   return(1);
}

/******************************************************************************
* Function: write_pgm_image
* Purpose: This function writes an image in PGM format. The file is either
* written to the file specified by outfilename or to standard output if
* outfilename = NULL. A comment can be written to the header if coment != NULL.
******************************************************************************/
int write_pgm_image(char *outfilename, unsigned char *image, int rows,
    int cols, char *comment, int maxval)
{
   FILE *fp;

   /***************************************************************************
   * Open the output image file for writing if a filename was given. If no
   * filename was provided, set fp to write to standard output.
   ***************************************************************************/
   if(outfilename == 0) fp = stdout;
   else{
      if((fp = fopen(outfilename, "w")) == 0){
         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
            outfilename);
         return(0);
      }
   }

   /***************************************************************************
   * Write the header information to the PGM file.
   ***************************************************************************/
   fprintf(fp, "P5\n%d %d\n", cols, rows);
   if(comment != 0)
      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
   fprintf(fp, "%d\n", maxval);

   /***************************************************************************
   * Write the image data to the file.
   ***************************************************************************/
   if((unsigned int)rows != fwrite(image, cols, rows, fp)){
      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
      if(fp != stdout) fclose(fp);
      return(0);
   }

   if(fp != stdout) fclose(fp);
   return(1);
}
