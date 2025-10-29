/****************************************************************/
/* Implementation of Ray Tracing Algorithm by Robert L. Siddon	*/
/* Refer to Medical Physics, 12(2), Mar/April 1985 pp.252-255	*/
/* pixel structure is not important, we can pretty much do      */
/* everything with just ray structure				*/
/* the pixel structure was taken out on July 30, 1993		*/
/****************************************************************/
/* catch the smallest and largest alpha's for the blurring */
/* July 11, 1994 */
/*****************/

#include <stdio.h>
#include <math.h>
// #include "math2.h"
#include <stdlib.h>
#include <string.h>
#include <malloc.h>
// #include "nrutil.h"
#include <string.h>

// USED FOR INTERFACE WITH MATLAB!!!!!!
#include "mex.h"
#define MAXNUMELEM 14     // max number of elements for the function (matlab)
#define BUFFPARALENG 30   // length of the buffer used to store input (in char)
#define NAMEFUNCTION "weigen3d"  // name of the function
/********************************/
/* spect definition header file */
/********************************/
#define	NPXL	64	/* # of pixels in each direction */
#define NDET	64	/* # of detectors in lateral and axil directions */
#define NPLANE	NDET	/* # of planes in projection */
#define NANGLE	64	/* # of projection angles, has to be EVEN number */
#define NRAYS	NDET*NANGLE /* number of total rays */
#define	PXL_HWD	0.712	/* 0.712 cm to 64x64 object; 45.568=0.712*64 */
#define	RAY_HWD	0.712	/* 0.356 cm to 128x128 object; 45.568=0.356*128 */
#define TRUE	1
#define FALSE	0
#define ROR	28.5	/* 29 cm, max error of 0.000190735 cm in weights */
                        /* 25 cm, max error of 0.529838    cm in weights */
                        /* 28 cm,  max error of 0.343189    cm in weights */
			/* 30 cm, max error of 0.000202179 cm in weights */
#define MAGIC	0.6999  /* trying to minimize the -1 plane blurring */
#define min(A,B)		((A)>(B)? (B) : (A))
#define max(A,B)		((A)>(B)? (A) : (B))
#define default_dir "./"

/* weights structure */
typedef struct wgt_stru { /* 16 bytes */
  float wgt;			/* 4 bytes */
  unsigned short int ipxl;	/* 2 bytes */
  unsigned short int vray;	/* 2 bytes */
  struct wgt_stru *next;/* 8 bytes */
} *wgt_ptr;

/* ray structure */
typedef struct ray_stru { /* 24 bytes */
  int pixels;	/* 4 bytes, number of pixels intercepting the ray */
  float dcp;    /* 4 bytes, distance from camera to roi */
  wgt_ptr hd_ptr; /* 8 bytes, points to the header of the weights list */ 
} *ray_ptr;

typedef struct alpha_stru { /* structure is defined for merging x, y's alpha */
  double alpha;		    /* 8 bytes */
  struct alpha_stru *next;  /* 8 bytes */
} *alpha_ptr; /* 16 bytes in total */

// #include "default.h"

void weigen3d(argc,argv)
  int argc;
  char *argv [MAXNUMELEM];
  
{
  int npxl,ndet,nangle,nplane,nangle_2,nangle_4,nrays,nrays_2;
  float r_major,r_minor,ror; /* radius */
  double *x,*y,*w; /* x, y and w are x, y and projection coordinates */
  double *ct,*st; /* vectors holding sin & cos of proj. angles */
  float pxl_wd,ray_wd,plane_wd;
  double x_l,x_r,y_d,y_u; /* x-leftmost, x-rightmost, y-downmost and y-upmost */
  int i,j,k;
  double theta,dtheta;
  double x1,y1,x2,y2; /* end pts of a ray in Siddon's */
  double ctr_x1,ctr_y1,ctr_x2,ctr_y2;/* coord. of center pt. of projection */
  unsigned char vertical,horizontal; /* flags of vertical/horizontal a ray */
  double ax_1,ax_nx,ay_1,ay_ny; /* alpha_x_1, _x_nx, _y_1 and _y_ny */
  double amin,amax; /* alpha_min and alpha_max */
  double abmin,abmax; /* alpha_min and alpha_max right at the object boundary */
  int imin,imax,jmin,jmax;
  double ax,ay;
  alpha_ptr a_header,a_p,a_q;
  alpha_ptr add_alpha();
  alpha_ptr sort_alpha();
  /* allocate an int vector with subscript range v[nl..nh] */
  unsigned short int **roi;
  /* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
  ray_ptr ray_vector();
  ray_ptr ray_header;
  int ix,iy,pix,piy;
  double xm,ym;
  double d12,dist,pdist;
  char fn[80],fn_roi[80],strbuf[80],strbuf1[80];
  FILE *fp;
  int ray_num;
  wgt_ptr wp,wq;
  double p_alpha,q_alpha,pp_alpha,pq_alpha;
  int same_pxl;
  float f_one=1,f_zero=0;
  float dcp,acc; /* Distance between Collimator to first Pixel */

  /* ----------------------------------- */
  /* command parsing and initializations */
  /* ----------------------------------- */
  /* printf("sizeof(alpha_ptr)= %-2d\n",sizeof(*a_header)); */
  /* printf("sizeof(ray_ptr)= %-2d\n",sizeof(*ray_header)); */
  /* printf("sizeof(wgt_ptr)= %-2d\n",sizeof(*wp)); */
  if(argc == 1 ) {
    printf("\nweigen3d for proj3d.\n");
    printf("data structures for 3d weigen are different from");
    printf(" the ones in 2d weigenp.\n");
printf("specify -np planes to sub-divide 2-D pixels for planar projection.\n");
    printf("  default has number of planes same as number of detectors.\n");
    printf("  if plane is set to 1 (np=1), weigen3d will be the same as ");
    printf("weigen2d.\n\n");
  printf("Usage: %s -p pxls[%-d] -d dets[%-d] -s steps[%-d] -np planes[%-d]\n",
         argv[0],NPXL,NDET,NANGLE,NPLANE);
    printf("       -rpx [reconstruction pixel size, %5.3fcm]\n",PXL_HWD);
    printf("       -apx [=rpx, acquisition pixel size, %5.3fcm]\n",RAY_HWD);
    printf("       [-roi field of view constraint for reconstruction space]\n");
    /* printf("       -ror [%-2dcm, radius of rotation]\n",ROR); */
    printf("\n");
    //exit(1);
	return;
  }
  npxl=NPXL;     /* initial number of pixels */
  ndet=NDET;     /* initial number of detectors */
  nplane=NPLANE; /* initual number of planes */
  nangle=NANGLE; /* initial number of angles */
  pxl_wd=PXL_HWD;
  ror=ray_wd=0;
  strcpy(fn_roi,"NULL");
  for(i=1; i<argc; i++) {
    if (strcmp(argv[i],"-p")==0)	/* number of pixels */
      sscanf(argv[i+1],"%d",&npxl), printf("npxl=%-d\n",npxl);
    else if (strcmp(argv[i],"-d")==0)	/* number of ndets */
      sscanf(argv[i+1],"%d",&ndet), printf("ndet=%-d\n",ndet);
    else if (strcmp(argv[i],"-s")==0)	/* number of steps */
      sscanf(argv[i+1],"%d",&nangle);
    else if (strcmp(argv[i],"-np")==0)	/* number of planes */
      sscanf(argv[i+1],"%d",&nplane);
    else if (strcmp(argv[i],"-roi")==0) /* region of interest */
      sscanf(argv[i+1],"%s",fn_roi);
    else if (strcmp(argv[i],"-rpx")==0)	/* reconstruction pixel size */
      sscanf(argv[i+1],"%f",&pxl_wd), printf("pxl_wd=%f\n",pxl_wd);
    else if (strcmp(argv[i],"-apx")==0)	/* acquisition pixel size */
      sscanf(argv[i+1],"%f",&ray_wd), printf("ray_wd=%f\n",ray_wd);
    else if (strcmp(argv[i],"-ror")==0)/* major radius */
      sscanf(argv[i+1],"%f",&ror);
  }
  if (ray_wd == 0) ray_wd=pxl_wd;
  plane_wd= (ndet/nplane)*ray_wd;
  if (ror == 0) ror=ROR;
  r_major=(npxl/2-MAGIC)*pxl_wd; /* note: ror has nothing to do with r_major */
  r_minor=r_major;
  nrays=ndet*nangle;
  nrays_2=nrays/2;
  nangle_2=nangle/2;
  nangle_4=nangle/4;
printf("#pixels= %-d, #detectors= %-d, #steps= %-d, #rays= %-d, #planes= %-d\n",
	        	npxl,ndet,nangle,nrays,nplane);
  printf(" pxl_wd = %9.7fcm, det_wd = %9.7fcm, plane_wd = %9.7fcm\n",
						    pxl_wd,ray_wd,plane_wd);
  printf(" rec. diameter = %6.3fcm", pxl_wd*npxl);
  printf(" r_major = %6.3fcm, r_minor = %6.3fcm\n",r_major,r_minor);
  printf(" radius of rotation is %6.3fcm\n",ror);
  if (strcmp(fn_roi,"NULL")==0) printf("default large ROI will be used\n");
  else printf("user defined ROI (%s) will be used\n",fn_roi);

  x =(double *) malloc(npxl*sizeof(double));
  y =(double *) malloc(npxl*sizeof(double));
  w =(double *) malloc(ndet*sizeof(double));
  ct=(double *) malloc(nangle*sizeof(double));
  st=(double *) malloc(nangle*sizeof(double));

  /* ===initialization=== */
  for (i=0; i<npxl; i++) { /* coordinates for 1-D pixels */
    x[i]=y[i]=(i+0.5-npxl/2)*pxl_wd; /* max 22.428 and 22.606 */
  }
  for (i=0; i<ndet; i++) { /* coordinates for 1-D projections */
    w[i]=(i+0.5-ndet/2)*ray_wd;
  }
  for (i=0,theta=0.0,dtheta=2*M_PI/nangle; i<nangle; i++) {
    theta=i*dtheta;
    ct[i]=cos(theta);
    st[i]=sin(theta);
  }

  /* index starting from 0 at lower left corner */
  /* 4032-> 4095 */
  /*    :->    : */
  /*   64->  127 */
  /*    0->   63 */
  roi= (unsigned short int **) usimatrix(0,npxl-1,0,npxl-1);
  if (strcmp(fn_roi,"NULL") == 0) {
    for (j=0; j<npxl; j++) {
      for (i=0; i<npxl; i++) {
        if (((x[i]*x[i])/(r_major*r_major)+(y[j]*y[j])/(r_minor*r_minor))<=1.0)
         roi[j][i]= 1; /* then it is in field of view */
        else roi[j][i] = 0;
      }
    }
  }
  else {
	  //has to be modified!!!!!! ROI is defined!!!!!
    strcpy(fn,fn_roi);
    fp=fopen(fn,"rb");
    for (j=0; j<npxl; j++)
      fread(roi[npxl-j-1],npxl*sizeof(unsigned short int),1,fp);
    fclose(fp);
  }

  x_l=x[0]-pxl_wd/2.0; /* locates the left -most x plane coordinate */
  x_r=x[npxl-1]+pxl_wd/2.0; /* locates the right-most x plane coordinate */
  y_d=y[0]-pxl_wd/2.0; /* locates the down -most y plane coordinate */
  y_u=y[npxl-1]+pxl_wd/2.0; /* locates the up   -most y plane coordinate */
  printf("x_l=%5.2f, x_r=%5.2f, y_d=%5.2f, y_u=%5.2f\n",x_l,x_r,y_d,y_u);

  /* ================================ */
  /* prepare a string of ????.???x??? */
  strcpy(strbuf,"_");
  strcat(strbuf,"p");
  sprintf(strbuf1,"%-d",npxl);
  strcat(strbuf,strbuf1);
  strcat(strbuf,"d");
  sprintf(strbuf1,"%-d",ndet);
  strcat(strbuf,strbuf1);
  strcat(strbuf,"s");
  sprintf(strbuf1,"%-d",nangle);
  strcat(strbuf,strbuf1);
  strcat(strbuf,"np");
  sprintf(strbuf1,"%-d",nplane);
  strcat(strbuf,strbuf1);
  strcat(strbuf,"_3D");
  /* ================================ */
  strcpy(fn,default_dir);
  strcat(fn,"rec"); /* stands for record */
  strcat(fn,strbuf);
  fp=fopen(fn,"w");
  fprintf(fp,"#pixels = %3d, #detectors = %3d, #steps = %3d, #rays = %5d,",
		npxl,ndet,nangle,nrays);
  fprintf(fp," #planes= %-d\n",nplane);
  fprintf(fp," pxl_wd = %9.7fcm, det_wd = %9.7fcm, plane_wd = %9.7fcm\n", 
		pxl_wd,ray_wd,plane_wd);
  fprintf(fp," rec. diameter = %6.3fcm, r_major = %6.3fcm, r_minor = %6.3fcm\n",
  	        pxl_wd*npxl,r_major,r_minor);
  /* fprintf(fp," radius of rotation = %6.3fcm\n",ror); */
  fclose(fp);
  /* ================= */
  strcpy(fn,default_dir);
  strcat(fn,".rec"); /* stands for record */
  strcat(fn,strbuf);
  fp=fopen(fn,"wb");
  fwrite(&pxl_wd,sizeof(float),1,fp);
  fwrite(&ray_wd,sizeof(float),1,fp);
  fwrite(&ror,sizeof(float),1,fp);
  fwrite(&plane_wd,sizeof(float),1,fp);
  fclose(fp);
  /* ================= */
  /* ================================ */
  strcpy(fn,default_dir);
  strcat(fn,"roi"); /* stands for ray driven */
  strcat(fn,strbuf);
  fp=fopen(fn,"wb");
  for (j=0; j<npxl; j++)
    fwrite(roi[npxl-j-1],npxl*sizeof(unsigned short int),1,fp);
  fclose(fp);

  ray_header= ray_vector(0,nrays_2-1); /* <-------important structure */

  /* --------------------------------------------------------- */
  /* SPECT acquisition starts at 0, which is -PI/2 in geometry */
  /* --------------------------------------------------------- */
  for (i=ray_num=0,theta=0.0,dtheta=2*M_PI/nangle,d12=ror*2; i<nangle_2; i++) {
    theta=i*dtheta;
    if ((i == 0) || (i == nangle_2)) vertical=TRUE;
    				else vertical=FALSE;
    if ((i == nangle_4) || (i == nangle_4+nangle_2)) horizontal=TRUE;
						else horizontal=FALSE;
    ctr_x2=ror*cos(theta-M_PI_2); 
    ctr_y2=ror*sin(theta-M_PI_2);
    ctr_x1=ror*cos(theta+M_PI_2); 
    ctr_y1=ror*sin(theta+M_PI_2);

    for (j=0; j<ndet; j++,ray_num++) {

      /* ================== */
      /* ray tracing starts */
      /* ================== */
      x2=ctr_x2+w[j]*ct[i]; /* X2 of (2) in Siddon's */
      y2=ctr_y2+w[j]*st[i]; /* Y2 of (2) in Siddon's */
      x1=ctr_x1+w[ndet-j-1]*ct[(i+nangle_2)%nangle]; /* X1 of (2) in Siddon's */
      y1=ctr_y1+w[ndet-j-1]*st[(i+nangle_2)%nangle]; /* Y1 of (2) in Siddon's */
      /* d12=sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1)); */

      ray_header[ray_num].pixels= 0;
      ray_header[ray_num].hd_ptr= NULL;

      ax_1=ax_nx=ay_1=ay_ny=0.0;
      if (!vertical) /* i.e. x2 != x1 */ {
	if (x2 > x1) {
	  ax_1 =(x_l-x1)/(x2-x1); /* (4) in Siddon's */
	  ax_nx=(x_r-x1)/(x2-x1);
	  if (ax_1 > ax_nx) printf("woo! fatal error! ax_1 > ax_nx");
	}
	else /* x2 < x1 */ {
	  ax_1 =(x_r-x1)/(x2-x1); /* (4) in Siddon's (added by tsp) */
	  ax_nx=(x_l-x1)/(x2-x1);
	  if (ax_1 > ax_nx) printf("woo! fatal error! ax_1 > ax_nx");
	}
      }
      if (!horizontal) /* i.e. y2 !=y1 */ {
	if (y2 > y1) {
	  ay_1 =(y_d-y1)/(y2-y1); /* (4) in Siddon's */
	  ay_ny=(y_u-y1)/(y2-y1);
	  if (ay_1 > ay_ny) printf("woo! fatal error! ay_1 > ay_ny");
	}
	else /* y2 < y1 */ {
	  ay_1 =(y_u-y1)/(y2-y1); /* (4) in Siddon's */
	  ay_ny=(y_d-y1)/(y2-y1);
	  if (ay_1 > ay_ny) printf("woo! fatal error! ay_1 > ay_ny");
	}
      }
      /* ================================================= */
      /* to implement equation (5) in Siddon's		   */
      /* if [vertical] then ax_1 and ax_nx are undefined   */ 
      /* if [horizontal] then ay_1 and ay_ny are undefined */
      /* no ray should be horizontal and vertical	   */
      /* ================================================= */
      if ((!horizontal) && (!vertical)) {
	amin=max(max(0.0,ay_1),ax_1);
      	amax=min(min(1.0,ay_ny),ax_nx);
      }
      else if (horizontal) {
	if ((y2 < y_d) || (y2 > y_u)) {
	  amin=1.0;
	  amax=0.0;
	}
	else {
	  amin=max(0.0,ax_1);
      	  amax=min(1.0,ax_nx);
	}
      }
      else /* must be vertical */ {
	if ((x2 < x_l) || (x2 > x_r)) {
	  amin=1.0;
	  amax=0.0;
	}
	else {
	  amin=max(0.0,ay_1);
      	  amax=min(1.0,ay_ny);
	}
      }

      /* ================================================= */
      /* to implement equation (6) in Siddon's		   */
      /* if [vertical] then ax_1 and ax_nx are undefined   */ 
      /* if [horizontal] then ay_1 and ay_ny are undefined */
      /* no ray should be horizontal and vertical	   */
      /* planes are numbered from 0 to 64(npxl)		   */
      /* ================================================= */
      if (amin < amax) { /* has some interception */
	a_header = NULL;
	if (!vertical) {
	  if (x2 > x1) {
	    /* ============================================================== */
	    /* DON'T do double-to-integer conversion directly 		      */
	    /* DO double-to-float first, THEN float-to-integer as shown below */
	    /* double-to-integer may give [7 to 57], not [8 to 57]            */
	    /* ============================================================== */
	    imin=(float) (0.5+(x1+amin*(x2-x1))/pxl_wd+npxl/2);
	    imax=(float) (0.5+(x1+amax*(x2-x1))/pxl_wd+npxl/2);
	    /* printf("imin=%4d, imax=%4d\n",imin,imax); */
	    for (k=imax; k>=imin; k--) {/* descending order */
	      ax=((k-npxl/2)*pxl_wd-x1)/(x2-x1); 
		/* ^^^^^^^^^^^^^^^^ locates the plane coordinate */
	      a_header=add_alpha(a_header,ax);
				  /* make it ascending order in list */ 
	    }
	  }
	  else {
	    imin=(float) (0.5+(x1+amax*(x2-x1))/pxl_wd+npxl/2);
	    imax=(float) (0.5+(x1+amin*(x2-x1))/pxl_wd+npxl/2);
	    /* printf("imin=%4d, imax=%4d\n",imin,imax); */
	    for (k=imin; k<=imax; k++) { /* descending order */
	      ax=((k-npxl/2)*pxl_wd-x1)/(x2-x1); 
		/* ^^^^^^^^^^^^^^^^ locates the plane coordinate */
	      a_header=add_alpha(a_header,ax);
				  /* make it ascending order in list */ 
	    }
	  }
	}
	else imin=imax=(-1);

	if (!horizontal) {
	  if (y2 > y1) {
	    jmin=(float) (0.5+(y1+amin*(y2-y1))/pxl_wd+npxl/2);
	    jmax=(float) (0.5+(y1+amax*(y2-y1))/pxl_wd+npxl/2);
	    /* printf("jmin=%4d, jmax=%4d\n",jmin,jmax); */
	    for (k=jmax; k>=jmin; k--) {
	      ay=((k-npxl/2)*pxl_wd-y1)/(y2-y1); 
		/* ^^^^^^^^^^^^^^^^ locates the plane coordinate */
	      if (vertical) a_header=add_alpha(a_header,ay);
					  /* a_header still NULL */
	      else a_header=sort_alpha(a_header,ay);
	    }
	  }
	  else {
	    jmin=(float) (0.5+(y1+amax*(y2-y1))/pxl_wd+npxl/2);
	    jmax=(float) (0.5+(y1+amin*(y2-y1))/pxl_wd+npxl/2);
	    /* printf("jmin=%4d, jmax=%4d\n",jmin,jmax); */
	    for (k=jmin; k<=jmax; k++) {
	      ay=((k-npxl/2)*pxl_wd-y1)/(y2-y1); 
		/* ^^^^^^^^^^^^^^^^ locates the plane coordinate */
	      if (vertical) a_header=add_alpha(a_header,ay);
					  /* a_header still NULL */
	      else a_header=sort_alpha(a_header,ay);
	    }
	  }
	}
	else jmin=jmax=(-1);

	/* circular FOV projection */
	for (a_p=a_header,a_q=a_p->next,abmin=1,abmax=0;
	     a_q!=NULL; a_p=a_p->next,a_q=a_q->next) {
	  if ((a_q->alpha -a_p->alpha) < 0.0)
	  {
			printf("Severe Error\n");
			return;
	  }
	  /* (xm,ym) represents the middle point of two intersections */
	  xm= x1+(a_p->alpha*(x2-x1)+a_q->alpha*(x2-x1))/2.0;
	  ym= y1+(a_p->alpha*(y2-y1)+a_q->alpha*(y2-y1))/2.0;
          if ((xm >= x_l) && (xm <= x_r) && (ym >= y_d) && (ym <= y_u)) {
	    ix=(float) floor((xm/pxl_wd+npxl/2)); /* added on 10-2-95 */
	    iy=(float) floor((ym/pxl_wd+npxl/2)); /* added on 10-2-95 */
	    dist=d12*(a_q->alpha-a_p->alpha);
	    p_alpha=a_p->alpha;
	    q_alpha=a_q->alpha;
	    same_pxl=FALSE;
	    if (a_p != a_header) {
	      if ((ix==pix) && (iy==piy)) {
		same_pxl=TRUE;
		/*
		printf("same_pxl =%8.5f p=%8.5f\n",dist,pdist);
		printf("p1= %8.5f, q1= %8.5f, p2 =%8.5f, q2= %8.5f\n",
			p_alpha,q_alpha,pp_alpha,pq_alpha);
		*/
	      }
	    }
	    pix=ix;
	    piy=iy;
	    pdist=dist;
	    pp_alpha=p_alpha;
	    pq_alpha=q_alpha;
	    /* printf("(%7.4f, %7.4f)->(%4d, %4d)\n",xm,ym,ix,iy); */
	    if (same_pxl) {
	      abmin=min(abmin,p_alpha);
	      abmax=max(abmax,q_alpha);
	      wq->wgt += (float) d12*(a_q->alpha - a_p->alpha);
	    }
	    else if (roi[iy][ix] == 1) {
	      /* ========================================================= */
	      /* start to prepare backward (projection to pixel) At matrix */
	      /* ========================================================= */
	      abmin=min(abmin,p_alpha);
	      abmax=max(abmax,q_alpha);
	      ++(ray_header[ray_num].pixels);
	      wq=(wgt_ptr)malloc(sizeof(*wq));
	      wq->wgt = (float) d12*(a_q->alpha - a_p->alpha);
	      wq->ipxl = iy*npxl+ix; /* one dimentional image array */
	      wq->vray = -1;
	      wq->next = NULL;
	      if (ray_header[ray_num].hd_ptr == NULL) {
		ray_header[ray_num].hd_ptr= wp = wq;
	      }
	      else { /* close to far */
		/* wp->next = wq;
		   wp=wq;  */
		wq->next = ray_header[ray_num].hd_ptr;
		ray_header[ray_num].hd_ptr = wq;
	      }
	    }
          }
	}
	/* circular FOV projection */
	
	ray_header[ray_num].dcp=(1-abmax)*d12; /* close */
	/* ray_header[ray_num].last_d=abmin*d12; */ 
        /* free the list for alpha parameters */
	for (a_p=a_header; a_p!=NULL; a_q=a_p,a_p=a_p->next,free(a_q));
      }
    } /* for j loop */
  } /* for i loop */

  /***************************************************/
  /* to store At matrix (backward projection matrix) */
  /***************************************************/
  strcpy(fn,default_dir);
  strcat(fn,"wgt");
  strcat(fn,strbuf);
  fp=fopen(fn,"wb");
  printf("storing At (backward) matrix: %s...\n",fn);
  for (i=0; i<nrays_2; i++) {
    if (ray_header[i].pixels == 0) {
      fwrite(&(ray_header[i].pixels),sizeof(int),1,fp);
      fwrite(&(f_zero),sizeof(float),1,fp);
    }
    else { /* major revision starts here */
      wp=ray_header[i].hd_ptr;
      dcp=ray_header[i].dcp-(ror-(ndet/2)*ray_wd);
      if (dcp < 0) { /* examine possible cause */
	printf("pixel out bnd i.e. dcp(%7.5f)<0 and reset dcp to 0...\n",dcp);
	dcp=0;
      }

      for (j=(float) (dcp/plane_wd),acc=dcp-plane_wd*j;
	  (j<nplane)&&(wp!=NULL); wp=wp->next) {
	if (((wp->wgt)+acc) > plane_wd) { /* split the weight to two pixels */
	  wq=(wgt_ptr)malloc(sizeof(*wq));
	  ++ray_header[i].pixels;
	  (*wq) = (*wp); 
	  wp->next = wq;
	  wp->wgt = plane_wd-acc;
	  wq->wgt -= wp->wgt;
	  wp->vray = j;
	  acc=0;
	  ++j; /* index to the next plane */
	}
	else { /* no split, and stay in the same plane */
	  acc += wp->wgt; 
	  wp->vray=j;
	}
      }
      fwrite(&(ray_header[i].pixels),sizeof(int),1,fp);
      fwrite(&(dcp),sizeof(float),1,fp);
      /*
      if (i<(3*ndet)) printf(" now has [%-d] pixles and dcp=%7.2f",
      ray_header[i].pixels,dcp);
      */
      for (j=0, wp=ray_header[i].hd_ptr; wp != NULL;
	wq=wp, wp=wp->next, free(wq), j++) {
	/* the above statement was modified to alleviate the memory demand
	   during the weigen3d */
	fwrite(&(wp->ipxl),sizeof(unsigned short int),1,fp);
	fwrite(&(wp->vray),sizeof(unsigned short int),1,fp);
	fwrite(&(wp->wgt),sizeof(float),1,fp);
	if (wp->vray == (-1)) printf("ray#[%-d] vray#[%-d]\n",i,wp->vray);
      }
      if (j !=ray_header[i].pixels)
	nrerror("Error at storing matrix A!");
    }
  }
  fclose(fp);
  printf(" done!\n");

} /* main */

/* ----------------------------- */
/* put alpha into list structure */
/* ----------------------------- */
alpha_ptr add_alpha(header,alpha)
  alpha_ptr header;
  double alpha;
{
  alpha_ptr p;
  
  p= (alpha_ptr) malloc(sizeof(*header)); /* allocate space for the pixel */
  p->alpha= alpha;
  p->next= NULL;
  if (header == NULL) header=p;
  else p->next= header;

  return(p);
}

/* ------------------------------ */
/* sort alpha into list structure */
/* ------------------------------ */
alpha_ptr sort_alpha(header,alpha)
  alpha_ptr header;
  double alpha;
{
  alpha_ptr p,q,r;

  r= (alpha_ptr) malloc(sizeof(*header)); /* allocate space for the pixel */
  r->alpha= alpha;
  r->next= NULL;
  if (header == NULL) return(r);
  else {
    for (p=NULL, q=header; (q != NULL) && (alpha > q->alpha); p=q, q=q->next); 
    if (p==NULL) {
      r->next=header;
      return(r);
    }
    else {
      p->next=r;
      if (q != NULL) r->next=q;
      return(header);
    }
  }
}

/* ------------------------------------------ */
/* create a linear array of ray_stru elements */
/* ------------------------------------------ */
ray_ptr ray_vector(nl,nh)
int nl,nh;
{
        ray_ptr v;

        v=(ray_ptr) malloc((unsigned) (nh-nl+1)*sizeof(*v));
        if (!v) nrerror("allocation failure in ray_vector()");
        return v-nl;
}
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@ MATLAB FUNCTION INTERFACE  @@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])

{
	int j,buflen;		// number of arguments
	double *input1;		// buffer for input numbers
	double intpar;		// used to store integer part, not used....
	double remai;		// used to store fractional portion
	char *argvp [MAXNUMELEM+1];	// input table given to C function
	char buff[BUFFPARALENG];	// buffer
	
	// check even number of arguments
	if ((nrhs%2)!=0 && nrhs!=1)
		mexErrMsgTxt("You don't have entered an even number of parameters\n Try to be more careful next time!!");
	if (nrhs>14)
		mexErrMsgTxt("Too many arguments!!!\n The maximum number of arguments allowed is 14!\n Try to be more careful next time!!");
	// fill the table
	if (nrhs!=1)
	{
		for (j=0;j<(nrhs/2);j++)
		{
			// FIRST PART OF THE EXTRACTION
			// extraction of the string
			// check if string
			if (mxIsChar(prhs[2*j])!=1)
				mexErrMsgTxt("You messed up the arguments!!!\nThe odd numbered arguments must be STRINGS!!!!\n Try to be more careful next time!!");
			// get length of input string
			buflen=(mxGetM(prhs[2*j])*mxGetN(prhs[2*j]))+1;
			// allocation of the memory
			argvp[2*j+1]=mxCalloc(buflen,sizeof(char));
			//put string
			mxGetString(prhs[2*j],argvp[2*j+1],buflen);
		
			// SECOND PART OF THE EXTRACTION
			//check if double
			if ((mxIsDouble(prhs[2*j+1]))==1)
			{
				// get the arguments
				input1=mxGetPr(prhs[2*j+1]);
				remai=modf(*input1,&intpar);
				if (remai==0)
				{
					// case integer
					sprintf(buff,"%d",(int) *input1);
					// check length and allocate memory
					buflen=strlen(buff)+1;
					argvp[2*j+2]=mxCalloc(buflen,sizeof(char));
					// put string in input table
					strcpy(argvp[2*j+2],buff);
				}
				else
				{
					//case float
					sprintf(buff,"%f",(float) *input1);
					// check length and allocate memory
					buflen=strlen(buff)+1;
					argvp[2*j+2]=mxCalloc(buflen,sizeof(char));
					// put string in input table
					strcpy(argvp[2*j+2],buff);
				}
			}
			else
			{
				// case string
				// get length of string
				buflen=(mxGetM(prhs[2*j+1])*mxGetN(prhs[2*j+1]))+1;
				// allocate memory
				argvp[2*j+2]=mxCalloc(buflen,sizeof(char));
				//put string
				mxGetString(prhs[2*j+1],argvp[2*j+2],buflen);
			}	
		}
	}

	// fill unutilized elements of input table by NULL vector
	for (j=nrhs;j<MAXNUMELEM;j++)
	{
		argvp[j+1]=NULL;
	}
	// fill first element with null
	// in C the first element of argv is the name of the function
	argvp[0]=NULL;
	
	if (nrhs==1)
	{
		buflen=strlen(NAMEFUNCTION)+1;
		argvp[0]=mxCalloc(buflen,sizeof(char));  // the name of the function has 8+'\0' char
		strcpy(argvp[0],NAMEFUNCTION);   // gives name of the function, will be displayed in the help
		argvp[1]=NULL;
	}
	// call C function
	weigen3d(nrhs,argvp);

}
