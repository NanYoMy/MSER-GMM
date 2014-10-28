#include<stdlib.h>
#include<string.h>
#include<assert.h>
#include<fstream>
#include <vector>
#include <queue>
#include <set>
#include <iostream>
using namespace std;

typedef unsigned int idx_t ;
typedef  unsigned short val_t;
typedef  int vl_bool;
#define VL_MSER_VOID_NODE 0xffffffff
#define MAX_PIXEL_VALUE 65535
#define VL_MSER_PIX_MAXVAL 65536
#define VL_MAX(x,y) (((x)>(y))?(x):(y))
#define VL_MIN(x,y) (((x)<(y))?(x):(y))
#define TEST
#define ROW 16

#define COL 16

class Point{

public:
	float col;
	float row;
	float value;
	Point(){
		
		col=0;
		row=0;
		value=-1;

	}

};



struct _VlMserReg
{
	idx_t parent ;   /**< points to the parent region.            */


	idx_t shortcut ; /**< points to a region closer to a root. ���ָ��һ��extremal region   */
	idx_t height ;   /**< region height in the forest.            */
	idx_t area ;     /**< area of the region.                     */
} ;

/** @internal @brief MSER: basic region */
typedef struct _VlMserReg VlMserReg ;

/* ----------------------------------------------------------------- */
/** @internal @brief MSER: extremal region (declaration)
**
** Extremal regions (ER) are extracted from the region forest. Each
** region is represented by an instance of this structure. The
** structures are stored into an array, in arbitrary order.
**
** ER are arranged into a tree. @a parent points to the parent ER, or
** to itself if the ER is the root.
**
** An instance of the structure represents the extremal region of the
** level set of intensity VlMserExtrReg::value and containing the
** pixel VlMserExtReg::index.
**
** VlMserExtrReg::area is the area of the extremal region and
** VlMserExtrReg::area_top is the area of the extremal region
** containing this region in the level set of intensity
** VlMserExtrReg::area + @c delta.
**
** VlMserExtrReg::variation is the relative area variation @c
** (area_top-area)/area.
**
** VlMserExtrReg::max_stable is a flag signaling whether this extremal
** region is also maximally stable.
**/
struct _VlMserExtrReg
{
	int          parent ;     /**< index of the parent region                   */
	int          index ;      /**< index of pivot pixel                         */
	val_t		value ;      /**< value of pivot pixel                         */
	idx_t      shortcut ;   /**< shortcut used when building a tree           */
	idx_t      area ;       /**< area of the region                           */
	float    variation ;  /**< rel. area variation                          */
	idx_t      max_stable ; /**< max stable number (=0 if not maxstable)      */

	int maxBranch;
	std::set<int>* MSERchildrens;
	bool MSERNonOverlap ;

	_VlMserExtrReg(){

		MSERNonOverlap=true;
		maxBranch=0;
		MSERchildrens=new std::set<int>;
	}

	~_VlMserExtrReg(){
	
		delete MSERchildrens;

	}


} ;
typedef struct _VlMserExtrReg VlMserExtrReg ;

struct _VlMserStats
{
	int num_extremal ;      /**< number of extremal regions                                */
	int num_unstable ;      /**< number of unstable extremal regions                       */
	int num_abs_unstable ;  /**< number of regions that failed the absolute stability test */
	int num_too_big ;       /**< number of regions that failed the maximum size test       */
	int num_too_small ;     /**< number of regions that failed the minimum size test       */
	int num_duplicates ;    /**< number of regions that failed the duplicate test          */
} ;

typedef struct _VlMserStats VlMserStats ;

struct _VlMserFilt
{

	/** @name Image data and meta data @internal */
	/*@{*/
	int                ndims ;   /**< number of dimensions                    */
	int               *dims ;    /**< dimensions                              */
	int                nel ;     /**< number of image elements (pixels)       */
	int               *subs ;    /**< N-dimensional subscript                 */
	int               *dsubs ;   /**< another subscript                       */
	int               *strides ; /**< strides to move in image data           */
	/*@}*/

	idx_t           *perm ;    /**< pixel ordering                          */
	idx_t           *joins ;   /**< sequence of join ops                    */
	int                njoins ;  /**< number of join ops                      */

	/** @name Regions */
	/*@{*/
	VlMserReg         *r ;       /**< basic regions                           */
	VlMserExtrReg     *er ;      /**< extremal tree                           */
	idx_t           *mer ;     /**< maximally stable extremal regions       */
	int                ner ;     /**< number of extremal regions              */
	int                nmer ;    /**< number of maximally stable extr. reg.   */
	int                rer ;     /**< size of er buffer                       */
	int                rmer ;    /**< size of mer buffer                      */
	/*@}*/

	/** @name Ellipsoids fitting */
	/*@{*/
	double         *acc ;     /**< moment accumulator.                    */
	double         *ell ;     /**< ellipsoids list.                       */
	int                rell ;    /**< size of ell buffer                     */
	int                nell ;    /**< number of ellipsoids extracted         */
	int                dof ;     /**< number of dof of ellipsoids.           */

	/*@}*/

	/** @name Configuration */
	/*@{*/
	vl_bool   verbose ;          /**< be verbose                             */
	int       delta ;            /**< delta filter parameter                  */
	double    max_area ;         /**< badness test parameter                 */
	double    min_area ;         /**< badness test parameter                 */
	double    max_variation ;    /**< badness test parameter                 */
	double    min_diversity ;    /**< minimum diversity                      */
	/*@}*/

	VlMserStats stats ;          /** run statistic                           */
} ;

typedef struct _VlMserFilt VlMserFilt ;

inline bool isOneBranch(const VlMserExtrReg* er, const int i){

	bool flag=false;
	int index=i;

	while (1==er[index].MSERchildrens->size())
	{
		int tmp=*(er[index].MSERchildrens->begin());
		index=tmp;
	}

	if (er[index].MSERchildrens->size()==0)
	{
		return true;
	}else{
		return false;
	}
}

inline Point getPoint(idx_t index){

	idx_t row=ROW;
	idx_t col=COL;
	Point p;
	int temp=index/col;
	p.row=temp;
	p.col = index-(temp)*col;
	return p;
}




void
outPutExreamalRigons(const val_t *image,std::vector<idx_t>& er )
{


	/************************************************************************/
	/* ͨ�����������vector��er����Ӧ��region��������                                                                    */
	/************************************************************************/
	std::queue<idx_t> q;
	std::vector<std::vector<idx_t>> extreamRegions;
	std::vector<idx_t> er_pt=er;
	int k, nel, ndims ; 
	int  * dims;
	val_t const * I_pt =image;
	int last = 0 ;
	int last_expanded = 0 ;
	val_t value = 0 ;
	int*   subs_pt ;       /* N-dimensional subscript                 */
	int*   nsubs_pt ;      /* diff-subscript to point to neigh.       */
	idx_t* strides_pt ;    /* strides to move in image array          */
	idx_t* visited_pt ;    /* flag                                    */

	/* get dimensions */
	ndims = 2;
	dims=(int*)malloc(sizeof(int)*ndims);
	dims[0]=ROW;
	dims[1]=COL;
	nel   =ROW*COL;
	


	/* allocate stuff */
	subs_pt    =(int*) malloc( sizeof(int)      * ndims ) ;
	nsubs_pt   =(int*) malloc( sizeof(int)      * ndims ) ;
	strides_pt =(idx_t*) malloc( sizeof(idx_t)    * ndims ) ;
	visited_pt =(idx_t*) malloc( sizeof(idx_t)    * nel   ) ;


	/* compute strides to move into the N-dimensional image array */
	strides_pt [0] = 1 ;
	for(k = 1 ; k < ndims ; ++k) {
		strides_pt [k] = strides_pt [k-1] * dims [k-1] ;
	}

	/* load first pixel */
	idx_t erSize=er_pt.size();
	for(idx_t i=0; i<erSize ;i++ ){

		memset(visited_pt, 0, sizeof(idx_t) * nel) ;
		std::vector<idx_t> region;
		//initial the queue
		idx_t rootIndex=er_pt.at(i)-1;
		q.push(rootIndex);
		region.push_back(rootIndex);
		value=I_pt[rootIndex];
		visited_pt [rootIndex] = 1 ;
		//���������
		while(!q.empty()){


			/* pop next node xi */
			idx_t index = q.front();
			q.pop();

			/* convert index into a subscript sub; also initialize nsubs 
			to (-1,-1,...,-1) */
			{
				idx_t temp = index ;
				for(k = ndims-1 ; k >=0 ; --k) {
					nsubs_pt [k] = -1 ;
					subs_pt  [k] = temp / strides_pt [k] ;
					temp         = temp % strides_pt [k] ;
				}
			}

			/* process 8 neighbors of xi */
			while( true ) {
				int good = true ;
				idx_t nindex = 0 ;

				/* compute NSUBS+SUB, the correspoinding neighbor index NINDEX
				and check that the pixel is within image boundaries. */
				for(k = 0 ; k < ndims && good ; ++k) {
					int temp = nsubs_pt [k] + subs_pt [k] ;
					good &= 0 <= temp && temp < dims[k] ;
					nindex += temp * strides_pt [k] ;
				}      

				/* process neighbor
				1 - the pixel is within image boundaries;
				2 - the pixel is indeed different from the current node
				(this happens when nsub=(0,0,...,0));
				3 - the pixel has value not greather than val
				is a pixel older than xi
				4 - the pixel has not been visited yet
				*/
				if(good 
					&& nindex != index 
					&& I_pt [nindex] <= value
					&& ! visited_pt [nindex] ) {

						/* mark as visited */
						visited_pt [nindex] = 1 ;

						/* add to list */
						region.push_back(nindex);
						q.push(nindex);
				}

				/* move to next neighbor */      
				k = 0 ;
				while(++ nsubs_pt [k] > 1) {
					nsubs_pt [k++] = -1 ;
					if(k == ndims) goto done_all_neighbors ;
				}
			} /* next neighbor */
done_all_neighbors : ;

		}

		extreamRegions.push_back(region);
		region.clear();
	}

	/*
	* Save results
	*/
	{
		//all is saved in members_pt

		//std::ofstream otf("D:\\mser\\data\\result.txt");
		std::ofstream posAndIntensity("D:\\mser\\data\\posAndIntensity.txt");
		std::vector<Point> centers; 
		std::vector<Point> highestIntensity; 
		int extreamSize=extreamRegions.size();
		//posAndIntensity<<"*****************************"<<std::endl;
		for (int i=0;i<extreamSize;i++)
		{
			float totalIntensity=0;
			float colIntensity=0;
			float rowIntensity=0;
			float currentIntensity=0;
			Point p;
			std::vector<idx_t> tmp=extreamRegions.at(i);
			int tmpSize=tmp.size();
			Point hightIntensityPoint;
			
			for(int j=0;j<tmpSize;j++){

				p=getPoint(tmp.at(j));
				
				currentIntensity=MAX_PIXEL_VALUE-I_pt[tmp.at(j)];


				int tmpCol=static_cast<int>(p.col)+1;
				int tmpRow=static_cast<int>(p.row)+1;
				assert(tmpCol<=COL);

				/*��ȡ���ֵ*/
				if (currentIntensity>hightIntensityPoint.value)
				{
					hightIntensityPoint.value=currentIntensity;
					hightIntensityPoint.col=tmpCol;
					hightIntensityPoint.row=tmpRow;
				}else{}

			
				posAndIntensity<<tmpCol<<" "<<tmpRow<<" "<<currentIntensity<<std::endl;


			}
			
			highestIntensity.push_back(hightIntensityPoint);


		}
		posAndIntensity.close();

		ofstream out("D:\\mser\\data\\highestIntensity.txt");
		int highIntensitySize=highestIntensity.size();

		/************************************************************************/
		/* out put the center                                                    */
		/************************************************************************/

		for(int i=0;i<highIntensitySize;i++){

			out<<highestIntensity.at(i).col<<" "<<highestIntensity.at(i).row<<" "<<highestIntensity.at(i).value<<endl;
		}
		out.close();
	}


	

	/* free stuff */
	free( visited_pt ) ;
	free( strides_pt ) ;
	free( nsubs_pt   ) ;
	free( subs_pt    ) ;
}


void
outPutCandidateRigons(const val_t *image,std::vector<idx_t>& er )
{


	/************************************************************************/
	/* ͨ�����������vector��er����Ӧ��region��������                                                                    */
	/************************************************************************/
	std::queue<idx_t> q;
	std::vector<std::vector<idx_t>> extreamRegions;
	std::vector<idx_t> er_pt=er;
	std::set<idx_t> pixelSet;
	int k, nel, ndims ; 
	int  * dims;
	val_t const * I_pt =image;
	int last = 0 ;
	int last_expanded = 0 ;
	val_t value = 0 ;
	int*   subs_pt ;       /* N-dimensional subscript                 */
	int*   nsubs_pt ;      /* diff-subscript to point to neigh.       */
	idx_t* strides_pt ;    /* strides to move in image array          */
	idx_t* visited_pt ;    /* flag                                    */

	/* get dimensions */
	ndims = 2;
	dims=(int*)malloc(sizeof(int)*ndims);
	dims[0]=ROW;
	dims[1]=COL;
	nel   =ROW*COL;

	/* allocate stuff */
	subs_pt    =(int*) malloc( sizeof(int)      * ndims ) ;
	nsubs_pt   =(int*) malloc( sizeof(int)      * ndims ) ;
	strides_pt =(idx_t*) malloc( sizeof(idx_t)    * ndims ) ;
	visited_pt =(idx_t*) malloc( sizeof(idx_t)    * nel   ) ;

	/* compute strides to move into the N-dimensional image array */
	strides_pt [0] = 1 ;
	for(k = 1 ; k < ndims ; ++k) {
		strides_pt [k] = strides_pt [k-1] * dims [k-1] ;
	}

	/* load first pixel */
	idx_t erSize=er_pt.size();
	for(idx_t i=0; i<erSize ;i++ ){

		memset(visited_pt, 0, sizeof(idx_t) * nel) ;
		std::vector<idx_t> region;
		//initial the queue
		idx_t rootIndex=er_pt.at(i)-1;
		q.push(rootIndex);
		region.push_back(rootIndex);
		value=I_pt[rootIndex];
		visited_pt [rootIndex] = 1 ;
		//���������
		while(!q.empty()){

			/* pop next node xi */
			idx_t index = q.front();
			q.pop();

			/* convert index into a subscript sub; also initialize nsubs 
			to (-1,-1,...,-1) */
			{
				idx_t temp = index ;
				for(k = ndims-1 ; k >=0 ; --k) {
					nsubs_pt [k] = -1 ;
					subs_pt  [k] = temp / strides_pt [k] ;
					temp         = temp % strides_pt [k] ;
				}
			}

			/* process 8 neighbors of xi */
			while( true ) {
				int good = true ;
				idx_t nindex = 0 ;

				/* compute NSUBS+SUB, the correspoinding neighbor index NINDEX
				and check that the pixel is within image boundaries. */
				for(k = 0 ; k < ndims && good ; ++k) {
					int temp = nsubs_pt [k] + subs_pt [k] ;
					good &= 0 <= temp && temp < dims[k] ;
					nindex += temp * strides_pt [k] ;
				}      

				/* process neighbor
				1 - the pixel is within image boundaries;
				2 - the pixel is indeed different from the current node
				(this happens when nsub=(0,0,...,0));
				3 - the pixel has value not greather than val
				is a pixel older than xi
				4 - the pixel has not been visited yet
				*/
				if(good 
					&& nindex != index 
					&& I_pt [nindex] <= value
					&& ! visited_pt [nindex] ) {

						/* mark as visited */
						visited_pt [nindex] = 1 ;

						/* add to list */
						region.push_back(nindex);
						q.push(nindex);
				}

				/* move to next neighbor */      
				k = 0 ;
				while(++ nsubs_pt [k] > 1) {
					nsubs_pt [k++] = -1 ;
					if(k == ndims) goto done_all_neighbors ;
				}
			} /* next neighbor */
done_all_neighbors : ;

		}
		extreamRegions.push_back(region);
		region.clear();
	}

	/*
	* Save results
	*/
	{

		int extreamSize=extreamRegions.size();
		//posAndIntensity<<"*****************************"<<std::endl;
		for (int i=0;i<extreamSize;i++)
		{
			float totalIntensity=0;
			float colIntensity=0;
			float rowIntensity=0;
			float currentIntensity=0;
			Point p;
			std::vector<idx_t> tmp=extreamRegions.at(i);
			int tmpSize=tmp.size();

			for(int j=0;j<tmpSize;j++){

					pixelSet.insert(tmp.at(j));
			}
		}
	}
	
	float currentIntensity=0;
	Point p;
	std::ofstream candidateRegions("D:\\mser\\data\\largetMSERRegion.txt");
	
	set<idx_t>::iterator endItr=pixelSet.end();

	for (set<idx_t>::iterator itr=pixelSet.begin();itr!=endItr;itr++)
	{
		p=getPoint(*itr);
		currentIntensity=MAX_PIXEL_VALUE-I_pt[*itr];
		
		int tmpCol=static_cast<int>(p.col)+1;
		int tmpRow=static_cast<int>(p.row)+1;
		assert(tmpCol<=COL);
		candidateRegions<<tmpCol<<" "<<tmpRow<<" "<<currentIntensity<<std::endl;

	}
	candidateRegions.close();

	/* free stuff */
	free( visited_pt ) ;
	free( strides_pt ) ;
	free( nsubs_pt   ) ;
	free( subs_pt    ) ;
}




/** -------------------------------------------------------------------
** @brief Climb the region forest to reach a root
**
** The function climbs the regions forest @r starting from the node
** @idx to the corresponding root.
**
** To speed-up the operation, the function uses the
** VlMserReg::shortcut field to quickly jump to the root. After the
** root is reached, all the used shortcut are updated.
**
** @param r regions' forest.
** @param idx stating node.
** @return index of the reached root.
**/

idx_t climb (VlMserReg* r, idx_t idx)
{

	idx_t prev_idx = idx ;
	idx_t next_idx ;
	idx_t root_idx ;

	/* move towards root to find it */
	/**/
	while (1) {

		/* next jump to the root */
		next_idx = r [idx] .shortcut ;

		/* recycle shortcut to remember how we came here */
		r [idx] .shortcut = prev_idx ;

		/* stop if the root is found */
		if( next_idx == idx ) break ;

		/* next guy */
		prev_idx = idx ;
		idx      = next_idx ;
	}

	root_idx = idx ;

	/*reset the short cut*/
	/* move backward to update shortcuts */
	while (1) {

		/* get previously visited one */
		prev_idx = r [idx] .shortcut ;

		/* update shortcut to point to the new root */
		r [idx] .shortcut = root_idx ;

		/* stop if the first visited node is reached */
		if( prev_idx == idx ) break ;

		/* next guy */
		idx = prev_idx ;
	}

	return root_idx ;
}

/** -------------------------------------------------------------------
** @brief Create a new MSER filter
**
** Initializes a new MSER filter for images of the specified
** dimensions. Images are @a ndims -dimensional arrays of dimensions
** @a dims.
**
** @param ndims number of dimensions.
** @param dims  dimensions.
**/

VlMserFilt* vl_mser_new (int ndims, int const* dims)
{
	VlMserFilt* f ;
	int *strides, k ;

	f =(VlMserFilt*) malloc (sizeof(VlMserFilt)) ;

	f-> ndims   = ndims ;
	f-> dims    = (int*)malloc (sizeof(int) * ndims) ;
	/*2��int��*/
	f-> subs    = (int*)malloc (sizeof(int) * ndims) ;
	f-> dsubs   = (int*)malloc (sizeof(int) * ndims) ;
	f-> strides = (int*)malloc (sizeof(int) * ndims) ;

	/* shortcuts */
	strides = f-> strides ;

	/* copy dims to f->dims */
	for(k = 0 ; k < ndims ; ++k) {
		//����dims 1024 1024
		f-> dims [k] = dims [k] ;
	}

	/* compute strides to move into the N-dimensional image array */
	strides [0] = 1 ;
	for(k = 1 ; k < ndims ; ++k) {

		strides [k] = strides [k-1] * dims [k-1] ;
	}

	//strides[0]=1 strides[1]=1024

	/* total number of pixels */
	//stides[k]=strides[k-1]*dims[k-1]
	f-> nel = strides [ndims-1] * dims [ndims-1] ;

	/* dof of ellipsoids */
	f-> dof = ndims * (ndims + 1) / 2 + ndims ;

	/* more buffers */
	/*����Ԫ�ظ����Ĵ�С*/
	f-> perm   = (idx_t*)malloc (sizeof(idx_t)   * f-> nel) ;
	f-> joins  = (idx_t*)malloc (sizeof(idx_t)   * f-> nel) ;
	f-> r      =(VlMserReg*) malloc (sizeof(VlMserReg) * f-> nel) ;

	f-> er     = 0 ;
	f-> rer    = 0 ;
	f-> mer    = 0 ;
	f-> rmer   = 0 ;
	f-> ell    = 0 ;
	f-> rell   = 0 ;

	/* other parameters */
	f-> delta         = 30 ;
	f-> max_area      = 0.75 ;
	f-> min_area      = 3.0 / f-> nel ;
	f-> max_variation = 0.25 ;
	f-> min_diversity = 0.5;

	return f ;
}

/** -------------------------------------------------------------------
** @brief Delete MSER filter
**
** The function releases the MSER filter @a f and all its resources.
**
** @param f MSER filter to be deleted.
**/

void vl_mser_delete (VlMserFilt* f)
{
	if(f) {
		if(f-> acc   )  ( f-> acc    ) ;
		if(f-> ell   )  ( f-> ell    ) ;

		if(f-> er    )  ( f-> er     ) ;
		if(f-> r     )  ( f-> r      ) ;
		if(f-> joins )  ( f-> joins  ) ;
		if(f-> perm  )  ( f-> perm   ) ;

		if(f-> strides) ( f-> strides) ;
		if(f-> dsubs  ) ( f-> dsubs  ) ;
		if(f-> subs   ) ( f-> subs   ) ;
		if(f-> dims   ) ( f-> dims   ) ;

		if(f-> mer    ) ( f-> mer    ) ;
		 (f) ;
	}
}


/** -------------------------------------------------------------------
** @brief Process image
**
** The functions calculates the Maximally Stable Extremal Regions
** (MSERs) of image @a im using the MSER filter @a f.
**
** The filter @a f must have been initialized to be compatible with
** the dimensions of @a im.
**
** @param f MSER filter.
** @param im image data.
**/

void vl_mser_process (VlMserFilt* f, val_t const* im)
{
	/* shortcuts */
	idx_t        nel     = f-> nel  ;
	idx_t       *perm    = f-> perm ;
	idx_t       *joins   = f-> joins ;
	int            ndims   = f-> ndims ;
	int           *dims    = f-> dims ;
	int           *subs    = f-> subs ;
	int           *dsubs   = f-> dsubs ;
	int           *strides = f-> strides ;
	VlMserReg     *r       = f-> r ;
	VlMserExtrReg *er      = f-> er ;
	idx_t       *mer     = f-> mer ;
	int            delta   = f-> delta ;

	int njoins = 0 ;
	int ner    = 0 ;
	int nmer   = 0 ;
	int nbig   = 0 ;
	int nsmall = 0 ;
	int nbad   = 0 ;
	int ndup   = 0 ;
	int overlap=0;

	int i, j, k ;

	/* delete any previosuly computed ellipsoid */
	f-> nell = 0 ;

	/* -----------------------------------------------------------------
	*   Sort pixels by intensity Ͱ����
	* -------------------------------------------------------------- */

	{
		idx_t buckets [ VL_MSER_PIX_MAXVAL ] ;

		/* clear buckets */
		memset (buckets, 0, sizeof(idx_t) * VL_MSER_PIX_MAXVAL ) ;

		/* compute bucket size (how many pixels for each intensity
		value) */
		for(i = 0 ; i < (int) nel ; ++i) {
			val_t v = im [i] ;
			++ buckets [v] ;
		}

		/* cumulatively add bucket sizes */
		for(i = 1 ; i < VL_MSER_PIX_MAXVAL ; ++i) {
			buckets [i] += buckets [i-1] ;
		}

		/* empty buckets computing pixel ordering */
		//��С��������perm[j]
		for(i = nel ; i >= 1 ; ) {

			//iԪ�ص�����ֵv
			val_t v = im [ --i ] ;
			//С������ֵv��Ԫ�صĸ���Ϊj����
			idx_t j = -- buckets [v] ;
			//i���Ԫ����perm��λ��Ϊj
			perm [j] = i ;
		}
	}

	/*set all the node's parent node to void*/
	/* initialize the forest with all void nodes */

	//���нڵ��parent��ֵΪVL_MSER_VOID_NODE
	for(i = 0 ; i < (int) nel ; ++i) {
		r [i] .parent = VL_MSER_VOID_NODE ;
	}


	/* -----------------------------------------------------------------
	*                        Compute regions and count extremal regions
	* -------------------------------------------------------------- */
	/*
	In the following:

	idx    : index of the current pixel
	val    : intensity of the current pixel
	r_idx  : index of the root of the current pixel
	n_idx  : index of the neighbors of the current pixel
	nr_idx : index of the root of the neighbor of the current pixel
	*/
	/* process each pixel by increasing intensity */
	for(i = 0 ; i < (int) nel ; ++i) {

		/* pop next node xi */
		//�õ�һ��Ԫ��
		idx_t     idx = perm [i] ;
		val_t val = im [idx] ;
		//��ǰ�ڵ�ĸ�
		idx_t     r_idx ;

		/* add the pixel to the forest as a root for now */
		r [idx] .parent   = idx ;
		r [idx] .shortcut = idx ;
		r [idx] .area     = 1 ;
		r [idx] .height   = 1 ;

		r_idx = idx ;

		/* convert the index IDX into the subscript SUBS; also initialize
		DSUBS to (-1,-1,...,-1) */
		{
			idx_t temp = idx ;
			for(k = ndims - 1 ; k >= 0 ; --k) {
				dsubs [k] = -1 ;
				subs  [k] = temp / strides [k] ;
				temp      = temp % strides [k] ;
			}
		}
		//subs[x][y]
		/* examine the neighbors of the current pixel */
		while (1) {
			idx_t n_idx = 0 ;
			vl_bool good = 1 ;

			/*
			Compute the neighbor subscript as NSUBS+SUB, the
			corresponding neighbor index NINDEX and check that the
			neighbor is within the image domain.
			*/
			for(k = 0 ; k < ndims && good ; ++k) {

				/*���Ͻǵ�neigbor*/
				int temp  = dsubs [k] + subs [k] ;
				good     &= (0 <= temp) && (temp < dims [k]) ;
				n_idx    += temp * strides [k] ;
			}

			/*
			The neighbor should be processed if the following conditions
			are met:

			1. The neighbor is within image boundaries.

			2. The neighbor is indeed different from the current node
			(the opposite happens when DSUB=(0,0,...,0)).

			3. The neighbor is already in the forest, meaning that it has
			already been processed.
			*/

			//r[n_idx]��r����region
			//parent�ڵ��������VL_MSER_VOID_NODE,����˵���ýڵ��Ѿ������ʡ�
			if (good &&
				n_idx != idx &&
				r [n_idx] .parent != VL_MSER_VOID_NODE ) {

					val_t	nr_val = 0 ;
					idx_t   nr_idx = 0 ;

					/*
					Now we join the two subtrees rooted at
					R_IDX = ROOT(  IDX)
					NR_IDX = ROOT(N_IDX).

					Note that R_IDX = ROOT(IDX) might change as we process more
					neighbors, so we need keep updating it.
					*/

					
					/*idx�ĸ������������е�Ԫ�ص�shortcutָ�����ĸ�*/
					r_idx = climb(r,   idx) ;
					//���ֵһ�����
					assert(im[r_idx]==val);

					/*idx�ĸ�*/
					nr_idx = climb(r, n_idx) ;

					int         hgt   = r [ r_idx] .height ;
					int         n_hgt = r [nr_idx] .height ;

					

					/*
					At this point we have three possibilities:

					��Ԫ���Ѿ���ͬһ��������
					(A) ROOT(IDX) == ROOT(NR_IDX). In this case the two trees
					have already been joined and we do not do anything.

					(B) I(ROOT(IDX)) == I(ROOT(NR_IDX)). In this case the pixel
					IDX is extending an extremal region with the same
					intensity value. Since ROOT(NR_IDX) will NOT be an
					extremal region of the full image, ROOT(IDX) can be
					safely added as children of ROOT(NR_IDX) if this
					reduces the height according to the union rank
					heuristic.

					(C) I(ROOT(IDX)) > I(ROOT(NR_IDX)). In this case the pixel
					IDX is starting a new extremal region. Thus ROOT(NR_IDX)
					WILL be an extremal region of the final image and the
					only possibility is to add ROOT(NR_IDX) as children of
					ROOT(IDX), which becomes parent.
					*/

					if( r_idx != nr_idx ) { /* skip if (A) */

						nr_val = im [nr_idx] ;

						if( nr_val == val && hgt < n_hgt ) {

							/* ROOT(IDX) becomes the child */
							/*��ȵ�ʱ��ֻҪ�������ص����Ϳ���*/
							r [r_idx]  .parent   = nr_idx ;
							//shortcut������Ľݾ�
							r [r_idx]  .shortcut = nr_idx ;
							r [nr_idx] .area    += r [r_idx] .area ;
							r [nr_idx] .height   = VL_MAX(n_hgt, hgt+1) ;

							joins [njoins++] = r_idx ;

						} else {

							/*���һ��ʼbuck�����ǰ�С���� ��ônr_val<=val���������buck�ǰ���С*/
							/* nr_val==val && hgt>n_hgt */
							/* nr_val<val */

							/* cases ROOT(IDX) becomes the parent */
							r [nr_idx] .parent   = r_idx ;
							r [nr_idx] .shortcut = r_idx ;
							r [r_idx]  .area    += r [nr_idx] .area ;
							r [r_idx]  .height   = VL_MAX(hgt, n_hgt + 1) ;

							/*���ϲ��ĵ�*/
							joins [njoins++] = nr_idx ;

							/*����ȵ�ʱ��region��+1 */
							/* count if extremal */
							if (nr_val != val) ++ ner ;

						} /* check b vs c */
					} /* check a vs b or c */
			} /* neighbor done */

			/* move to next neighbor */
			k = 0 ;
			/*����3*3��С��������������о���ǿ�������*/
			while(++ dsubs [k] > 1) {
				dsubs [k++] = -1 ;
				if(k == ndims) goto done_all_neighbors ;
			}
		} /* next neighbor */
done_all_neighbors:;
	} /* next pixel */

	/* the last root is extremal too */
	++ ner ;

	/* save back */
	f-> njoins = njoins ;

	/*��ֵ������*/
	f-> stats. num_extremal = ner ;

	/* -----------------------------------------------------------------
	*                                          Extract extremal regions
	* -------------------------------------------------------------- */

	/*
	Extremal regions are extracted and stored into the array ER.  The
	structure R is also updated so that .SHORTCUT indexes the
	corresponding extremal region if any (otherwise it is set to
	VOID).
	*/

	/* make room */
	if (f-> rer < ner) {
		if (er)  (er) ;
		//f->er  = er = (VlMserExtrReg*)malloc (sizeof(VlMserExtrReg) * ner) ;
		f->er=er=new VlMserExtrReg[ner];
		f->rer = ner ;
	} ;
	

	/* save back */
	f-> nmer = ner ;

	/* count again */
	ner = 0 ;

	/*��һ������extrel region!!!*/
	/* scan all regions Xi */
	for(i = 0 ; i < (int) nel ; ++i) {

		/* pop next node xi */
		idx_t     idx = perm [i] ;

		/*
		ÿһ������im[idx]����һ��r[idx]������Ӧ
		*/
		val_t val   = im [idx] ;
		idx_t     p_idx = r  [idx] .parent ;
		val_t p_val = im [p_idx] ;

		/* is extremal ? */
		/*���ڵ������ֵ���ڣ����ߵ����˸�*/
		vl_bool is_extr = (p_val > val) || idx == p_idx ;

		if( is_extr ) {

			/* if so, add it */
			/*extramal region ��index����Ϊ��ǰ�����ص��index*/
			er [ner] .index      = idx ;
	
	

			/*extreamal region�ĸ��ڵ��ȳ�ʼ�����Լ�*/
			er [ner] .parent     = ner ;
			er [ner] .value      = im [idx] ;
			er [ner] .area       = r  [idx] .area ;

			/* link this region to this extremal region */
			/*�ѵ�ǰ����ڵ���extrem region��ϵ����*/
			r [idx] .shortcut = ner ;

			/* increase count */
			++ ner ;
		} else {

			/* link this region to void */
			/*�Ǽ�ֵ�����ص㣬����һ�����ص�����ֵΪ A�� �Ѿ���һ�����ص������ֵΪA�ĳ�Ϊһ��extremal region*/
			r [idx] .shortcut =   VL_MSER_VOID_NODE ;
		}
	}

	/* -----------------------------------------------------------------
	*                                   Link extremal regions in a tree
	* -------------------------------------------------------------- */
	/*ÿ��er[i]��һ��extreamal region ��������ֵ����Ԫ�ص�*/

	/*һ���������ݽṹ��һ�����ɵ���������ɵ���*/
	/*һ������extremal region��ɵ���*/

	for(i = 0 ; i < ner ; ++i) {

		/*һ��er����һ�������������������er[i].index��ʾ�������ĸ�*/
		idx_t idx = er [i] .index ;

		/************************************************************************/
		/* Ѱ��                                                                     */
		/************************************************************************/
		do {
			idx = r[idx] .parent ;
			//���һ��r[k]��shortcutΪVL_MSER_VOID_NODE����ʾr[k]��parent����ֵ��r[k]һ��
		} while (r[idx] .shortcut == VL_MSER_VOID_NODE/*������shortcutΪVL_MSER_VOID_NODE��ʾΪ��*/) ;

		/*��������shortcut��2�����ã�������������ʱ��¼��ǰ�ĸ�*/
		/*������extream regionʱ��,��¼�������ڵ��Ӧ��extremal regions�ı��*/
		er[i] .parent   = r[idx] .shortcut ;

		/*extreamal region�����shortcut��ʾ��extream region id*/
		er[i] .shortcut = i ;
	}

	/* -----------------------------------------------------------------
	*                            Compute variability of +DELTA branches
	* -------------------------------------------------------------- */
	/* For each extremal region Xi of value VAL we look for the biggest
	* parent that has value not greater than VAL+DELTA. This is dubbed
	* `top parent'. */

	for(i = 0 ; i < ner ; ++i) {

		/* Xj is the current region the region and Xj are the parents */
		int     top_val = er [i] .value + delta ;
		int     top     = er [i] .shortcut ;

		/* examine all parents */
		/*�����Ż�����*/
		while (1) {
			int next     = er [top]  .parent ;
			int next_val = er [next] .value ;

			/* Break if:
			* - there is no node above the top or
			* - the next node is above the top value.
			*/
			if (next == top || next_val > top_val) break ;

			/* so next could be the top */
			top = next ;
		}

		/* calculate branch variation */
		{
			int area     = er [i  ] .area ;
			int area_top = er [top] .area ;
			er [i] .variation  = (float) (area_top - area) / area ;
			er [i] .max_stable = 1 ;
		}

		/* Optimization: since extremal regions are processed by
		* increasing intensity, all next extremal regions being processed
		* have value at least equal to the one of Xi. If any of them has
		* parent the parent of Xi (this comprises the parent itself), we
		* can safely skip most intermediate node along the branch and
		* skip directly to the top to start our search. */
		/************************************************************************/
		/* ����ط����Ż������⡣����                                  */
		/************************************************************************/
		{
			int parent = er [i] .parent ;
			int curr   = er [parent] .shortcut ;
			er [parent] .shortcut =  VL_MAX (top, curr) ;
		}
	}






	/* -----------------------------------------------------------------
	*                                  Select maximally stable branches
	* -------------------------------------------------------------- */
	printf("done\nExtremal regions: %d\n", ner) ;
	nmer = ner ;
	for(i = 0 ; i < ner ; ++i) {
		idx_t    parent = er [i     ] .parent ;
		val_t   val = er [i     ] .value ;
		float     var = er [i     ] .variation ;
		val_t p_val = er [parent] .value ;
		float   p_var = er [parent] .variation ;
		idx_t     loser ;

		/*
		Notice that R_parent = R_{l+1} only if p_val = val + 1. If not,
		this and the parent region coincide and there is nothing to do.
		*/

		/*����ط����ƺ����⡣Ӧ��p_val>=val+1�Ϳ��԰�*/
		if(p_val > val + 1) continue ;

		/* decide which one to keep and put that in loser */
		if(var < p_var) loser = parent ; else loser = i ;

		/* make loser NON maximally stable */
		if(er [loser] .max_stable) {
			-- nmer ;
			er [loser] .max_stable = 0 ;
		}
	}

#define SIXTEEN
#ifdef SIXTEEN
	std::ofstream outf("D:\\mser\\data\\16.txt");
	for (int i=0;i<ner;i++)
	{
		outf<<"ner: "<<i<<" value: "<<MAX_PIXEL_VALUE-(int)er[i].value<<" area"<<(int)er[i].area<<" parent"<<(int)er[i].parent<<std::endl;
		//outf<<"ner: "<<i<<" value: "<<(int)er[i].value<<" area"<<(int)er[i].area<<" parent"<<(int)er[i].parent<<std::endl;
	}

	outf.close();
#endif

	f-> stats. num_unstable = ner - nmer ;
	printf("mer region %d \n",nmer);

	/* -----------------------------------------------------------------
	*                                                 Further filtering
	* -------------------------------------------------------------- */
	/* It is critical for correct duplicate detection to remove regions
	* from the bottom (smallest one first).                          */
	{
		float max_area = (float) f-> max_area * nel ;
		float min_area = (float) f-> min_area * nel ;
		float max_var  = 1;//(float) f-> max_variation ;
		float min_div  = (float) f-> min_diversity ;

		/* scan all extremal regions (intensity value order) */
		for(i = ner-1 ; i >= 0L  ; --i) {

			/* process only maximally stable extremal regions */
			if (! er [i] .max_stable) continue ;
			if (er [i] .variation >= max_var ) { ++ nbad ;   goto remove ; }
			if (er [i] .area      >  max_area) { ++ nbig ;   goto remove ; }
			if (er [i] .area      <  min_area) { ++ nsmall ; goto remove ; }
			/*
			* Remove duplicates
			*/
			if (min_div < 1.0) {
				idx_t   parent = er [i] .parent ;
				int       area, p_area ;
				float div ;

				/* check all but the root mser */
				if((int) parent != i) {

					/* search for the maximally stable parent region */
					while(! er [parent] .max_stable) {
						idx_t next = er [parent] .parent ;
						if(next == parent) break ;
						parent = next ;
					}

					/* Compare with the parent region; if the current and parent
					* regions are too similar, keep only the parent. */
					area    = er [i]      .area ;
					p_area  = er [parent] .area ;
					div     = (float) (p_area - area) / (float) p_area ;

					if (div < min_div) { ++ ndup ; goto remove ; }
				} /* remove dups end */

			}
			continue ;
remove :
			er [i] .max_stable = 0 ;
			-- nmer ;
		} /* check next region */

		f-> stats .num_abs_unstable = nbad ;
		f-> stats .num_too_big      = nbig ;
		f-> stats .num_too_small    = nsmall ;
		f-> stats .num_duplicates   = ndup ;
	}


	/*���maxBranch construct the maximal stable extrem region tree*/
	for (int i=0;i<ner;i++)
	{

		idx_t up= i;
		int tmpChild=0;


		while (up!=er[up].parent)
		{

			/*�ҵ�һ��msere�ڵ�*/
			while ((1!=er[up].max_stable)&&(up!=er[up].parent))
			{
				up=er[up].parent;
			}

			/*�������ڵ�*/
			if (up==er[up].parent)
			{
				/*������max_stable*/
				if (1==er[up].max_stable)
				{
					continue;

				}else{}

				continue;
			}else{} //1==er[up].max_stable
			//
			tmpChild=up;
			up=er[up].parent;
			while((1!=er[up].max_stable)&&(up!=er[up].parent)){

				up=er[up].parent;
			}

			/*���ڵ�*/
			if (up==er[up].parent)
			{
				/*max_stable*/
				if (1==er[up].max_stable)
				{

#ifdef DEBUG
						std::cout<<"node: "<<up<<" get a child: "<<tmpChild<<std::endl;
#endif

				
					er[up].MSERchildrens->insert(tmpChild);

					er[up].maxBranch++;
				}else{ continue;}

			}else{
#ifdef DEBUG
				std::cout<<"node: "<<up<<" get a child: "<<tmpChild<<std::endl;
#endif

				er[up].MSERchildrens->insert(tmpChild);
				er[up].maxBranch++;
			}

			tmpChild=up;
			up=er[up].parent;
		}

	}


	printf("  Bad regions:        %d\n", nbad   ) ;
	printf("  Small regions:      %d\n", nsmall ) ;
	printf("  Big regions:        %d\n", nbig   ) ;
	printf("  Duplicated regions: %d\n", ndup   ) ;
	printf("  overlap regions: %d\n", overlap   ) ;
	printf("Cleaned-up regions: %d (%.1f%%)\n", 
		nmer, 100.0 * (double) nmer / ner) ;


#ifdef SIXTEEN
	std::ofstream outf2("D:\\mser\\data\\16-max.txt");
	for (int i=0;i<ner;i++)
	{
		if (er[i].max_stable)
		{
			outf2<<"ner: "<<i<<" value: "<<MAX_PIXEL_VALUE-(int)er[i].value<<" area"<<(int)er[i].area<<" parent"<<(int)er[i].parent<<std::endl;

		}
	}


	outf2.close();
#endif
	/* -----------------------------------------------------------------
	*                                                   Save the result
	* -------------------------------------------------------------- */

	/* make room */
	if (f-> rmer < nmer) {
		if (mer)  (mer) ;
		f->mer  = mer = (idx_t*)malloc( sizeof(idx_t) * nmer) ;
		f->rmer = nmer ;
	}

	/* save back */
	f-> nmer = nmer ;

	std::vector<idx_t> exIdx;

	j = 0 ;
	for (i = 0 ; i < ner ; ++i) {
		
		if ((er [i] .max_stable)&&(isOneBranch(er,i))){
			mer [j++] = er [i] .index ;
			exIdx.push_back(er[i].index+1);
		}
	}

	//out put the star regions 
	outPutExreamalRigons(im,exIdx);

	exIdx.clear();
	j = 0 ;
	for (i = 0 ; i < ner ; ++i) {

		if ((er [i] .max_stable)){
			mer [j++] = er [i] .index ;
			exIdx.push_back(er[i].index+1);
		}
	}

	//out put candidate regions
	outPutCandidateRigons(im,exIdx);

}


void main(){

	int* dims = (int*)malloc(sizeof(int)*2);
	dims[0]=ROW;
	dims[1]=COL;

	//initial the data
	val_t inPixel;
	val_t* image=new val_t[dims[0]*dims[1]];
	int a=0;
	std::ifstream inf("D:\\mser\\data\\candidate.txt");
	while(inf>>inPixel){

		image[a]=MAX_PIXEL_VALUE-inPixel;
		
		a++;
	}

	VlMserFilt* filter=vl_mser_new(2,dims);
	vl_mser_process(filter,image);
	vl_mser_delete(filter);

}