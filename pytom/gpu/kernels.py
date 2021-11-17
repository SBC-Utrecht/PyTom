spmvh_kernel_text = r"""

__device__ void atomic_add_float_twosum(float * address, float val, float * residue)
{
    // https://devtalk.nvidia.com/default/topic/817899/atomicadd-kahan-summation/
    // By Sylvain Collange
    // Also add the branch clsparse
    
    float x = atomicAdd(address, val);
    
    if (fabs(x) < fabs(val))
    {
        const float swap = x;
        x = val;
        val = swap;
    }
    float newacc = x + val;
    float r = val - (newacc - x);    // Recover rounding error using Fast2Sum, assumes |oldacc|>=|val|
    atomicAdd(residue, r);    // Accumulate error in residue
}    

__device__ void atomic_add_float( float *ptr, const float temp) 
{
    atomicAdd(ptr, temp); 
};   

__device__ void atomic_add_float2( float2 * ptr, const float2 temp, float2 *res) 
{
    atomic_add_float_twosum((float*)ptr, temp.x, (float*)res);
    atomic_add_float_twosum((float*)ptr+1, temp.y, (float*)res+1);
};  

extern "C" __global__
void spmvh(unsigned int Reps, unsigned int nRow, unsigned int prodJd, unsigned int sumJd, unsigned int dim, int *Jd, int *meshindex, int *kindx, 
           float2 *udata, float2 *k, float2 *res, float2 *input)
{

    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;
                                                                                                                                          
    unsigned int myRow0 = blockIdx.x*(blockSize) + tid;
    unsigned int myRow = myRow0/(float)Reps;
    unsigned int nc = myRow0 - myRow*Reps;
    float2 zero;
    zero.x = 0.0;
    zero.y = 0.0;
    //printf("%i %i %i\n", myRow, nRow, prodJd);
    if (myRow < nRow){ 
        for (unsigned int j = 0;  j  <  prodJd; j ++){
            float2 u = zero;

            // now doing the first dimension
            unsigned int index_shift = myRow * sumJd;
            //printf("rr: %i\n", index_shift);
            unsigned int J = Jd[0];
            unsigned int index =    index_shift +  meshindex[dim*j + 0];
            unsigned int col = kindx[index] ;
            float2 spdata = udata[index];
            index_shift += J; 
            //printf("rr2: %i\n", Jd[0]);
            for (unsigned int dimid = 1; dimid < dim; dimid ++ ){
                    J = Jd[dimid];
                    index =   index_shift + meshindex[dim*j + dimid];   // the index of the partial ELL arrays *kindx and *udata
                    col += kindx[index];// + 1  ;                                            // the column index of the current j
                    float tmp_x = spdata.x;
                    float2 tmp_udata = udata[index];
                    spdata.x = tmp_x * tmp_udata.x - spdata.y * tmp_udata.y;                            // the spdata of the current j
                    spdata.y = tmp_x * tmp_udata.y + spdata.y * tmp_udata.x; 
                    index_shift  += J;
            }; // Iterate over dimensions 1 -> Nd - 1
            //printf("rr3: %i\n", index_shift);

            float2 ydata = input[myRow*Reps + nc]; // kout[col];
            u.x =   spdata.x * ydata.x + spdata.y * ydata.y;
            u.y =  -spdata.y * ydata.x + spdata.x * ydata.y;
            //printf("tt: %i %f %f %i %i\n", myRow*Reps + nc, input[0].x, input[0].y, col, nc);
            atomic_add_float2((k + col*Reps + nc), u, (res + col*Reps + nc));

        }; // Iterate for (unsigned int j = 0;  j  <  prodJd; j ++)
    };  // if (m < nRow)
};    // End of pELL_spmvh_mCoil          
"""

spmv_kernel_text = r"""
extern "C" __global__ void spmv(    
        const    unsigned int   Reps,            // Number of coils
        const    unsigned int    nRow,        // number of rows
        const    unsigned int    prodJd,     // product of Jd
        const    unsigned int    sumJd,     // sum of Jd
        const    unsigned int    dim,           // dimensionality
        const unsigned int *Jd,            // Jd, length = dim
        const unsigned int *meshindex,            // meshindex, prodJd * dim
        const unsigned int *kindx,    // unmixed column indexes of all dimensions
        const float2 *udata,// interpolation data before Kronecker product
        float2 *vec,     // multi-channel kspace data, prodKd * Reps
        float2 *out)   // multi-channel output, nRow * Reps
{   
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;
    
    const unsigned int t = threadIdx.x; //blockIdx.x*(blockSize) + tid;
    const unsigned int vecWidth = 1024; //${LL}
    
    // Thread ID within wavefront
    const unsigned int id = t & (vecWidth-1);
    
    
    // One row per wavefront
    unsigned int vecsPerBlock = blockSize/vecWidth;
    unsigned int myRow = (blockIdx.x * vecsPerBlock) + (t/ vecWidth); // the myRow-th non-Cartesian sample
    unsigned int m = myRow / Reps;
    unsigned int nc = myRow - m * Reps;
    __shared__ float2 partialSums[1024]; // ${LL}
    float2 zero;
    zero.x = 0.0;
    zero.y = 0.0;
    partialSums[t] = zero;
    
    if (myRow < nRow * Reps)
    {
        const unsigned int vecStart = 0; 
        const unsigned int vecEnd = prodJd;             
        float2  y;//=zero;
        //printf("%i %i\n", id, vecEnd);
        for (unsigned int j = vecStart+id;  j < vecEnd; j += vecWidth)
        {    // now doing the first dimension
            unsigned int J = Jd[0];
            unsigned int index_shift = m * sumJd ;
            unsigned int index =    index_shift +  meshindex[dim*j + 0];
            unsigned int col = kindx[index] ;
            float2 spdata = udata[index];
            
            index_shift += J; 
            
            for (unsigned int dimid = 1; dimid < dim; dimid ++ )
            {
                unsigned int J = Jd[dimid];
                unsigned int index =  index_shift + meshindex[dim*j + dimid];   // the index of the partial ELL arrays *kindx and *udata
                col += kindx[index] ;//+ 1;                                            // the column index of the current j
                float tmp_x= spdata.x;
                float2 tmp_udata = udata[index];
                spdata.x = spdata.x * tmp_udata.x - spdata.y * tmp_udata.y;                            // the spdata of the current j
                spdata.y = tmp_x * tmp_udata.y + spdata.y * tmp_udata.x; 
                index_shift  += J;
            }
            //printf("%i ", col*Reps + nc);
            float2 vecdata = vec[col * Reps + nc];
            y.x =  spdata.x*vecdata.x - spdata.y*vecdata.y;
            y.y =  spdata.y*vecdata.x + spdata.x*vecdata.y;
            partialSums[t].x = y.x + partialSums[t].x;
            partialSums[t].y = y.y + partialSums[t].y;
            
        }
    
        __syncthreads(); 
    
        // Reduce partial sums
        unsigned int bar = vecWidth / 2;
        while(bar > 0)
        {
            if (id < bar)
            {
                partialSums[t].x = partialSums[t].x + partialSums[t+bar].x;
                partialSums[t].y = partialSums[t].y + partialSums[t+bar].y;
            }
            __syncthreads();
            
            bar = bar / 2;
        }            
        
        // Write result 
        if (id == 0)
        {
        out[myRow]=partialSums[t]; 
        }
    }
};  // End of pELL_spmv_mCoil
"""

sum_weighted_norm_complex_array_text = r'''
__device__ void warpReduceSum(volatile float* sdata, int tid, int blockSize) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
    if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
    if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
    if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
    if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



extern "C" __global__ 
void sum_weighted_norm_complex_array(float2 *g_idata, float *weight, float * g_odata, int n) {

    __shared__ float mean[1024];
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

    mean[tid] = 0.;                                                                                                                                                       

    while (i < n) {
        mean[tid] += weight[i] * (g_idata[i].x * g_idata[i].x + g_idata[i].y * g_idata[i].y); 
        if (i + blockSize < n){
                mean[tid] += weight[i+blockSize] * (g_idata[i + blockSize].x * g_idata[i + blockSize].x);
                mean[tid] += weight[i+blockSize] * (g_idata[i + blockSize].y * g_idata[i + blockSize].y);}
        i += gridSize;}
    __syncthreads();                                                                                                                                                       

    if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
    if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
    if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
    if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
    if (tid < 32){ warpReduceSum(mean, tid, blockSize);}                                                                                                                                                                      
    if (tid == 0) {g_odata[blockIdx.x] = mean[0];} 
    __syncthreads(); 

}'''

sum_norm_complex_array_text = r'''
__device__ void warpReduceSum(volatile float* sdata, int tid, int blockSize) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
    if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
    if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
    if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
    if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



extern "C" __global__ 
void sum_norm_complex_array(float2 *g_idata, float * g_odata, int n) {

    __shared__ float mean[1024];
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

    mean[tid] = 0.;                                                                                                                                                       

    while (i < n) {
        mean[tid] += (g_idata[i].x * g_idata[i].x + g_idata[i].y * g_idata[i].y); 
        if (i + blockSize < n){
                mean[tid] += (g_idata[i + blockSize].x * g_idata[i + blockSize].x);
                mean[tid] += (g_idata[i + blockSize].y * g_idata[i + blockSize].y);}
        i += gridSize;}
     __syncthreads();                                                                                                                                                       

    if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
    if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
    if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
    if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
    if (tid < 32){ warpReduceSum(mean, tid, blockSize);}                                                                                                                                                                      
    if (tid == 0) {g_odata[blockIdx.x] = mean[0];} 
   __syncthreads(); 

}'''

sum_text = r'''
__device__ void warpReduceSum(volatile float* sdata, int tid, int blockSize) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32];}
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16];}
    if (blockSize >= 16) {sdata[tid] += sdata[tid +  8];}
    if (blockSize >=  8) {sdata[tid] += sdata[tid +  4];}
    if (blockSize >=  4) {sdata[tid] += sdata[tid +  2];}
    if (blockSize >=  2) {sdata[tid] += sdata[tid +  1];} }



extern "C" __global__ 
void sumKernel(float *g_idata, float *g_mean,  int n) {

    __shared__ float mean[1024];
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

    mean[tid] = 0.;                                                                                                                                                       

    while (i < n) {
        mean[tid] += g_idata[i]; 
        if (i + blockSize < n){mean[tid] += g_idata[i + blockSize];}
        i += gridSize;}
     __syncthreads();                                                                                                                                                       

    if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512];} __syncthreads(); }                                                                                                    
    if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256];} __syncthreads(); }                                                                          
    if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128];} __syncthreads(); }                                                                          
    if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64];} __syncthreads(); }                                                                          
    if (tid < 32){ warpReduceSum(mean, tid, blockSize);}                                                                                                                                                                      
    if (tid == 0) {g_mean[blockIdx.x] = mean[0];} 
   __syncthreads(); 

}'''

cTensorCopy_text = '''
extern "C" __global__
void cTensorCopy(const unsigned int batch, const unsigned int dim, const  unsigned int *Nd_elements,
                            const  unsigned int *Kd_elements, const  float *invNd, const float2 *indata,
                            float2 *outdata, const int direction) {  
    
    
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    const unsigned int gid = blockIdx.x*(blockSize) + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2; 
    
    unsigned int curr_res = gid;
    unsigned int new_idx = 0;
    unsigned int group;
    
    for (unsigned int dimid =0; dimid < dim; dimid ++){
        group = (float)curr_res*invNd[dimid];
        new_idx += group * Kd_elements[dimid];
        curr_res = curr_res - group * Nd_elements[dimid];
    };
    
    if (direction == 1) {
        for (unsigned int bat=0; bat < batch; bat ++ )
        {
            outdata[new_idx*batch+bat]= indata[gid*batch+bat];
         };   
    };
    
    if (direction == -1) {
        for (unsigned int bat=0; bat < batch; bat ++ )
        {
            outdata[gid*batch+bat]= indata[new_idx*batch+bat];
        };   
    };
};'''

cTensorMultiply_text = """
extern "C" __global__
void cTensorMultiply(const unsigned int batch, const unsigned int dim, const  unsigned int *Nd, 
                                const unsigned int *Nd_elements, const  float *invNd_elements, const float *vec, 
                                float2 *outdata, const unsigned int div) 
{
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    const unsigned int gid = blockIdx.x*(blockSize) + tid; 
    const unsigned int pid = (float)gid / (float)batch;

    unsigned int group;
    unsigned int Nd_indx_shift = 0;
    float mul = 1.0; 
    unsigned int res = pid; 

    for (unsigned int dimid = 0; dimid < dim; dimid ++){
        group = (float)res * invNd_elements[dimid]; // The index along the axis
        res = res - group * Nd_elements[dimid];

        const unsigned int N = Nd[dimid]; 

        mul = mul * vec[group + Nd_indx_shift];

        Nd_indx_shift = Nd_indx_shift + N;
    }

    if (div == 1){
        // for (unsigned int bat = 0; bat < batch; bat ++ )
        // {
        float2 tmp = outdata[gid];
        tmp.x = tmp.x /  mul;
        tmp.y = tmp.y / mul;
        outdata[gid] = tmp;
        // };
    };
    
    if (div == 0){
        float2 tmp = outdata[gid];
        tmp.x = tmp.x *  mul;
        tmp.y = tmp.y * mul;
        outdata[gid] = tmp;
    };    

};"""

argmax_text = r'''

extern "C"  __device__ void warpReduce(volatile float* sdata, volatile int* maxid, int tid, int blockSize) {
    if (blockSize >= 64 && sdata[tid] < sdata[tid + 32]){ sdata[tid] = sdata[tid + 32]; maxid[tid] = maxid[tid+32];} 
    if (blockSize >= 32 && sdata[tid] < sdata[tid + 16]){ sdata[tid] = sdata[tid + 16]; maxid[tid] = maxid[tid+16];}
    if (blockSize >= 16 && sdata[tid] < sdata[tid +  8]){ sdata[tid] = sdata[tid +  8]; maxid[tid] = maxid[tid+ 8];}
    if (blockSize >=  8 && sdata[tid] < sdata[tid +  4]){ sdata[tid] = sdata[tid +  4]; maxid[tid] = maxid[tid+ 4];}
    if (blockSize >=  4 && sdata[tid] < sdata[tid +  2]){ sdata[tid] = sdata[tid +  2]; maxid[tid] = maxid[tid+ 2];}
    if (blockSize >=  2 && sdata[tid] < sdata[tid +  1]){ sdata[tid] = sdata[tid +  1]; maxid[tid] = maxid[tid+ 1];}
}

extern "C" __global__ 
void argmax(float *g_idata, float *g_odata, int *g_mdata, int n) {
    __shared__ float sdata[1024]; 
    __shared__ int maxid[1024];
    
    
    int blockSize = blockDim.x;                                                                                                                                                   
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

    sdata[tid] = -1000000;                                                          
    maxid[tid] = 0.;

    while (i < n) {
        //if (sdata[tid] < g_idata[i]){
        sdata[tid] = g_idata[i];
        maxid[tid] = i;
        //}
        if (i+blockSize < n){
            if (sdata[tid] < g_idata[i+blockSize]){
                sdata[tid] = g_idata[i+blockSize];
                maxid[tid] = i + blockSize;
            }
        }
        i += gridSize; 
    };
     __syncthreads();                                                                                                                                                       

    if (blockSize >= 1024){ if (tid < 512 && sdata[tid] < sdata[tid + 512]){ sdata[tid] = sdata[tid + 512]; maxid[tid] = maxid[tid+512]; } __syncthreads(); }
    if (blockSize >=  512){ if (tid < 256 && sdata[tid] < sdata[tid + 256]){ sdata[tid] = sdata[tid + 256]; maxid[tid] = maxid[tid+256]; } __syncthreads(); }
    if (blockSize >=  256){ if (tid < 128 && sdata[tid] < sdata[tid + 128]){ sdata[tid] = sdata[tid + 128]; maxid[tid] = maxid[tid+128]; } __syncthreads(); }
    if (blockSize >=  128){ if (tid <  64 && sdata[tid] < sdata[tid +  64]){ sdata[tid] = sdata[tid +  64]; maxid[tid] = maxid[tid+ 64]; } __syncthreads(); }
    if (tid < 32){ warpReduce(sdata, maxid, tid, blockSize);}                                                                                                                                                                      
    if (tid == 0) {g_odata[blockIdx.x] = sdata[0]; g_mdata[blockIdx.x] = maxid[0];} 


}'''

meanStdv_text = r'''
__device__ void warpReduce(volatile float* sdata, volatile float* stdv, int tid, int blockSize) {
    if (blockSize >= 64) {sdata[tid] += sdata[tid + 32]; stdv[tid] += stdv[tid+32];}
    if (blockSize >= 32) {sdata[tid] += sdata[tid + 16]; stdv[tid] += stdv[tid+16];}
    if (blockSize >= 16) {sdata[tid] += sdata[tid +  8]; stdv[tid] += stdv[tid+ 8];}
    if (blockSize >=  8) {sdata[tid] += sdata[tid +  4]; stdv[tid] += stdv[tid+ 4];}
    if (blockSize >=  4) {sdata[tid] += sdata[tid +  2]; stdv[tid] += stdv[tid+ 2];}
    if (blockSize >=  2) {sdata[tid] += sdata[tid +  1]; stdv[tid] += stdv[tid+ 1];} }

extern "C" __global__ 
void sumMeanStdv(float *g_idata, float *mask, float *g_mean, float * g_stdv,  int n) {
                                                                                                                                                     
    __shared__ float mean[1024];
    __shared__ float stdv[1024];
    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = blockIdx.x*(blockSize)*2 + tid;                                                                                                                       
    int gridSize = blockSize*gridDim.x*2;                                                                                                                         

    mean[tid] = 0.;                                                                                                                                                       
    stdv[tid] = 0.;
    
    while (i < n) {
        mean[tid] += g_idata[i] * mask[i] ; 
        stdv[tid] += g_idata[i] * g_idata[i] * mask[i];
        if (i + blockSize < n){
            mean[tid] += g_idata[i + blockSize] * mask[i+blockSize];
            stdv[tid] += g_idata[i + blockSize] * g_idata[i + blockSize] * mask[i+blockSize];}
        i += gridSize;}
     __syncthreads();                                                                                                                                                       

    if (blockSize >= 1024){ if (tid < 512) { mean[tid] += mean[tid + 512]; stdv[tid] += stdv[tid+512];} __syncthreads(); }                                                                                                    
    if (blockSize >= 512) { if (tid < 256) { mean[tid] += mean[tid + 256]; stdv[tid] += stdv[tid+256];} __syncthreads(); }                                                                          
    if (blockSize >= 256) { if (tid < 128) { mean[tid] += mean[tid + 128]; stdv[tid] += stdv[tid+128];} __syncthreads(); }                                                                          
    if (blockSize >= 128) { if (tid <  64) { mean[tid] += mean[tid +  64]; stdv[tid] += stdv[tid+ 64];} __syncthreads(); }                                                                          
    if (tid < 32){ warpReduce(mean, stdv, tid, blockSize);}                                                                                                                                                                      
    if (tid == 0) {g_mean[blockIdx.x] = mean[0]; g_stdv[blockIdx.x] = stdv[0];} 
   

}'''

reconstruction_wbp_text = '''

extern "C" __global__ 
void reconstruction_wbp(float * projection, float* proj_center, int * proj_dims, 
                        float * reconstruction, int * recon_center, int * recon_dims,
                        float * tr, int n) {

    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = (blockIdx.x*(blockSize) + tid);
    
    int x;
    int y;
    int z;
    
    //for(int d=0; d < recon_dims[2]; d++){
        if (i < n) {
            x = i  / (recon_dims[1] * recon_dims[2]) - recon_center[0]; 
            y = (i % (recon_dims[1] * recon_dims[2])) / recon_dims[2] - recon_center[1]; 
            z = (i % (recon_dims[1] * recon_dims[2])) % recon_dims[2] - recon_center[2];
            
            //float ix = tr[0] * float(x) + tr[2] * float(z) + float(proj_center[0]);
            //float iy = tr[4] * float(y) + float(proj_center[1]);
            
            float ix = tr[0] * float(x) + tr[1] * float(y) + tr[2] * float(z) + float(proj_center[0]);
			float iy = tr[3] * float(x) + tr[4] * float(y) + tr[5] * float(z) + float(proj_center[1]);
            
            int px = int(ix); 
            int py = int(iy);
    
            float xoff = ix - float(px); 
            float yoff = iy - float(py);
            
            if (px >= 1 && py >= 1 && px < proj_dims[0] && py < proj_dims[1]){
                float v1 = projection[(px-1)*proj_dims[1]+py-1] + (projection[px*proj_dims[1]+py-1]-projection[(px-1)*proj_dims[1]+py-1]) * xoff;
                float v2 = projection[(px-1)*proj_dims[1]+py  ] + (projection[px*proj_dims[1]+py  ]-projection[(px-1)*proj_dims[1]+py  ]) * xoff;
                float v3 = v1 + (v2-v1)*yoff;
                
                atomicAdd( reconstruction + i, v3);
                //atomicAdd( weights + i, 1);
            };
        };
        //i++;
    //};
};
'''

projection_text = '''

extern "C" __global__ 
void reconstruction_wbp_calc(float * projection, float* proj_center, int * proj_dims, 
                        float * reconstruction, int * recon_center, int * recon_dims,
                        float * tr, int m, int n) {

    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = (blockIdx.x*(blockSize) + tid);

    int x;
    int y;
    int z;

    //for(int d=0; d < recon_dims[2]; d++){
        if (i < n) {
            x = i  / (recon_dims[1] * recon_dims[2]) - recon_center[0]; 
            y = (i % (recon_dims[1] * recon_dims[2])) / recon_dims[2] - recon_center[1]; 
            z = (i % (recon_dims[1] * recon_dims[2])) % recon_dims[2] - recon_center[2];

            float ix = tr[0] * float(x) + tr[2] * float(z) + float(proj_center[0]) +0.5;
            float iy = tr[1] * float(y) + float(proj_center[1]) + 0.5;

            int px = int(ix); 
            int py = int(iy);

            float xoff = ix - float(px); 
            float yoff = iy - float(py);

            if (px >= 0 && px < proj_dims[0] && py > 0 && py < proj_dims[1]){
                atomicAdd( projection + px * proj_dims[1] + py, xoff*yoff*reconstruction[i]);

            if (px-1 >= 0 && px-1 < proj_dims[0] && py > 0 && py < proj_dims[1]){
                atomicAdd( projection + (px-1) * proj_dims[1] + py, (1-xoff)*yoff*reconstruction[i]);

            if (px-1 >= 0 && px-1 < proj_dims[0] && py-1 > 0 && py-1 < proj_dims[1]){
                atomicAdd( projection + (px-1)* proj_dims[1] + py-1, (1-xoff)*(1-yoff)*reconstruction[i]);

            if (px >= 0 && px < proj_dims[0] && py-1 > 0 && py-1 < proj_dims[1]){
                atomicAdd( projection + px * proj_dims[1] + py-1, xoff*(1-yoff)*reconstruction[i]);

            };
        };
        //i++;
    //};
};
'''

transformSpline_txt = '''

__device__ float interpolNearestNeighbour(float * data, float x, float y, float z, float defaultval, int sizex, int sizey, int sizez) {

    x += 0.5;
    y += 0.5;
    z += 0.5;
    if (x<0 || y<0 || z<0 || x>=sizex || y>=sizey || z>=sizez) {
        return defaultval;
    }else {
        return data[int(x) *  sizey* sizez + int(y) * sizez + int(z)];
    }
}

__device__ float interpolTriLinear(float *data, float x, float y, float z, float defaultval, int sizex, int sizey, int sizez) {

    if (x<0. || y<0. || z<0. || x > sizex || y > sizey || z > sizez ) {
        return defaultval;
    }
    int floorx = int(x);
    int floory = int(y);
    int floorz = int(z);
    float xoffseth = x - float(floorx);
    float yoffseth = y - float(floory);
    float zoffseth = z - float(floorz);

    int stridey = sizex;
    int stridez = sizex*sizey;

    float * src = &data[floorz * stridez + floory * stridey + floorx];

    switch (int(xoffseth > 0.) | (int(yoffseth > 0.)<<1) | (int(zoffseth > 0.)<<2)) {
        case 0x07: {
                float x00x, x01x, x10x, x11x, x0yx;
                x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += stridey;
                x01x = src[0] + (src[1]-src[0])*xoffseth;
                src += stridez;
                x11x = src[0] + (src[1]-src[0])*xoffseth;
                src -= stridey;
                x10x = src[0] + (src[1]-src[0])*xoffseth;
                x0yx = x00x + (x01x - x00x) * yoffseth;
                return x0yx + ((x10x + (x11x - x10x) * yoffseth) - x0yx) * zoffseth;
            }
        case 0x06: {
                float x0y0 = src[0] + (src[stridey]-src[0])*yoffseth;
                src += stridez;
                return x0y0 + (src[0] + (src[stridey]-src[0])*yoffseth - x0y0)*zoffseth;
            }
        case 0x05: {
                float x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += stridez;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)*zoffseth;
            }
        case 0x03: {
                float x00x = src[0] + (src[1]-src[0])*xoffseth;
                src += stridey;
                return x00x + (src[0] + (src[1]-src[0])*xoffseth - x00x)* yoffseth;
            }
        case 0x04:
            return src[0] + (src[stridez]-src[0])*zoffseth;
        case 0x02:
            return src[0] + (src[stridey]-src[0])*yoffseth;
        case 0x01:
            return src[0] + (src[1]-src[0])*xoffseth;
    }
    return src[0];
}



__device__ float CUB_INT(float x, float y, float xj, float x2, float x3){\
return y*(x-x2)/(xj-x2)*(x-x3)/(xj-x3);
}

__device__ float interpolTriCubic(float *data, float x, float y, float z, float defaultval, int sizex, int sizey, int sizez) {

    bool is3D = sizez > 1;
    int stridey = sizex;
    int stridez = sizex*sizey;


	// is the current position in the data or outside. return default value if outside
    if (is3D && (x<1. || y<1. || z<1. || x > (sizex)-2 || y > (sizey)-2 || z > (sizez)-2) ) {
    	return defaultval; //interpolTriLinear(data, x, y, z, defaultval, sizex, sizey, sizez);
    }

	if(!is3D && (x<1. || y<1. || x > (sizex)-2 || y > (sizey)-2) ){
		return defaultval; //interpolTriLinear(data, x, y, z, defaultval, sizex, sizey, sizez);
	}


    int floorx = (int(x));              //integer position along x
    int floory = (int(y));              //integer position along y
    int floorz = (int(z));              //integer position along z
    float xoffseth = x - float(floorx); //floating point offset along x
    float yoffseth = y - float(floory); //floating point offset along y
    float zoffseth = z - float(floorz); //floating point offset along z

    float * src = &data[int(floorz) * stridez + int(floory)*stridey + int(floorx)]; //point into the closest position to interpolation

    switch (int(xoffseth > 0.) | (int(yoffseth > 0.)<<1) | (int(zoffseth > 0.)<<2)) {
		case 0x00:{
                //all match grid, no interpolation
			return src[0];
		}
		default: {

			float px_pl_1 = floorx + 1; //all voxels plus 1 from x -> p3,p6,p9
			float px_mi_1 = floorx - 1; //all voxels minus 1 from x -> p1,p4,p7

			float py_pl_1 = floory + 1; //all voxels plus 1 from y -> p1,p2,p3
			float py_mi_1 = floory - 1; //all voxels minus 1 from x -> p7,p8,p9

			//those will hold the current voxel values
			float v1 = defaultval;
			float v2 = defaultval;
			float v3 = defaultval;
			float v4 = defaultval;
			float v5 = defaultval;
			float v6 = defaultval;
			float v7 = defaultval;
			float v8 = defaultval;
			float v9 = defaultval;

			float line1 = defaultval; //interpolation values for each line (y)
			float line2 = defaultval;
			float line3 = defaultval;

			int lowerLayerBound = -1;
			int upperLayerBound = 1;
			int layerOffset = 1;
			float layerValues[3];

			if(!is3D){
				lowerLayerBound = 0;
				upperLayerBound = 0;
				layerOffset = 0;
			}

			//interpolation values for each layer (z)
			for (int zIteration=lowerLayerBound; zIteration <= upperLayerBound; zIteration++) {
				//current position in memory plus current z layer offset in voxels (of type T)
				//first will be negative (-1), second 0 (same layer), third is 1, next layer

				//load the pixel values
				v1 = *(src-1-stridey + zIteration*stridez); //one line up in y direction, one position back in x
				v2 = *(src  -stridey + zIteration*stridez); //one line up in y direction
				v3 = *(src+1-stridey + zIteration*stridez); //one line up in y direction, one position forward in x
				v4 = *(src-1 + zIteration*stridez); //same line in y
				v5 = *(src + zIteration*stridez); //...
				v6 = *(src+1 + zIteration*stridez);
				v7 = *(src-1+stridey + zIteration*stridez);
				v8 = *(src  +stridey + zIteration*stridez);
				v9 = *(src+1+stridey + zIteration*stridez);

				//interpolate first row 1 2 3
				line1 = CUB_INT(x,v1,px_mi_1,floorx,px_pl_1); //px1,px2,px3
				line2 = CUB_INT(x,v2,floorx,px_pl_1,px_mi_1); //px2,px3,px1
				line3 = CUB_INT(x,v3,px_pl_1,px_mi_1,floorx); //px3,px1,px2
				//store values into v1
				v1 = line1+line2+line3;


				//same for the next rows
				line1 = CUB_INT(x,v4,px_mi_1,floorx,px_pl_1); //px4,px5,px6
				line2 = CUB_INT(x,v5,floorx,px_pl_1,px_mi_1); //px5,px6,px4
				line3 = CUB_INT(x,v6,px_pl_1,px_mi_1,floorx); //px6,px5,px4
				v2 = line1+line2+line3;

				line1 = CUB_INT(x,v7,px_mi_1,floorx,px_pl_1); //px7,px8,px9
				line2 = CUB_INT(x,v8,floorx,px_pl_1,px_mi_1); //px8,px9,px7
				line3 = CUB_INT(x,v9,px_pl_1,px_mi_1,floorx); //px9,px7,px8
				v3 = line1+line2+line3;

				//interpolate col 2 5 8 in y direction
				line1 = CUB_INT(y,v1,py_mi_1,floory,py_pl_1); //py2,py5,py8
				line2 = CUB_INT(y,v2,floory,py_pl_1,py_mi_1); //py5,py8,py2
				line3 = CUB_INT(y,v3,py_pl_1,py_mi_1,floory); //py8,py2,py5
				layerValues[zIteration + layerOffset] = line1+line2+line3;

			    //printf("GPU %d - %.1f %.1f %.1f - %.1f %.1f %.1f %.1f - %.1f %.1f %.1f %.1f - %.1f\\n", zIteration, x, y, z, v1, v2, v3, v4, line1, line2, line3, line4, layerValues[zIteration + layerOffset]);

			}

			if(is3D){
				line1 = CUB_INT(z,layerValues[0],(floorz-1),floorz,(floorz+1));
				line2 = CUB_INT(z,layerValues[1],floorz,(floorz+1),(floorz-1));
				line3 = CUB_INT(z,layerValues[2],(floorz+1),(floorz-1),floorz);
			}

			return (line1 + line2 + line3);

		}


		return defaultval;
    }
    return defaultval;
}

__device__ float CSPL_INT(float f1,float f2, float f3, float f4, float a, float b, float c, float d){
    return 3*(f1*(b-a)+f2*(c-a)+f3*(d-b)+f4*(d-c));
}

__device__ float CSPL_CALC(float c2, float c3, float D1, float D2, float off){
    return c2+D1*off+(3*(c3-c2)-2*D1-D2)*off*off+(2*(c2-c3)+D1+D2)*off*off*off;
}


__device__ float interpolCubicSpline(float * data, float x, float y, float z, float defaultval, int sizex, int sizey, int sizez) {

    int stridex = sizey;
    int stridey = sizex;
    int stridez = sizex*sizey;

    //printf("%d %d %d; \\n", sizex, sizey, sizez);

    bool is3D = sizez > 1;



	// is the current position in the data or outside. return default value if outside
    if (is3D && (x<2. || y<2. || z<2. || x > (sizex)-3 || y > (sizey)-3 || z > (sizez)-3) ) {
        return interpolTriCubic(data, x,y,z, defaultval, sizex, sizey, sizez);
    }

    if(!is3D && (x<2. || y<2. || x > (sizex)-3 || y > (sizey)-3) ){
        return interpolTriCubic(data, x,y,z, defaultval, sizex, sizey, sizez); 
    }

    int floorx = int(x);//integer position along x
    int floory = int(y);//integer position along y
    int floorz = int(z);//integer position along z
    float xoffseth = x - float(floorx);//floating point offset along x
    float yoffseth = y - float(floory);//floating point offset along y
    float zoffseth = z - float(floorz);//floating point offset along z

    float * src = &data[floorz * stridez + floory*stridey + floorx];//point into the closest position to interpolation

    switch (static_cast<int>(xoffseth > 0.) | (static_cast<int>(yoffseth > 0.)<<1) | (static_cast<int>(zoffseth > 0.)<<2)) {
		case 0x00:{
			//all match grid, no interpolation
			return src[0];
		}
		default:{
			/*interpolation square

			 p1   p2   p3 p4
			 p5   p6   p7 p8
			         P
			 p9   p10  p11 p12
			 p13  p14  p15 p16

			 */

			//stridey = stridex;

			float f1= -0.1556;
			float f2=  0.3111;
			float f3= -0.0889;
			float f4=  0.0444;

			//those will hold the current voxel values
			float v1 = defaultval;
			float v2 = defaultval;
			float v3 = defaultval;
			float v4 = defaultval;
			float v5 = defaultval;
			float v6 = defaultval;
			float v7 = defaultval;
			float v8 = defaultval;
			float v9 = defaultval;
			float v10 = defaultval;
			float v11 = defaultval;
			float v12 = defaultval;
			float v13 = defaultval;
			float v14 = defaultval;
			float v15 = defaultval;
			float v16 = defaultval;

			float line1 = defaultval; 
			float line2 = defaultval;
			float line3 = defaultval;
			float line4 = defaultval;

			float D1 = defaultval;
			float D2 = defaultval;

			int lowerLayerBound = -1;
			int upperLayerBound = 2;
			int layerOffset = 1;
			float layerValues[4];

			if(!is3D){
				lowerLayerBound = 0;
				upperLayerBound = 0;
				layerOffset = 0;
			}
			//interpolation values for each layer (z)
			for (int zIteration=lowerLayerBound; zIteration <= upperLayerBound; zIteration++) {

				//load the pixel values
				v1 = *(src-1-stridey + zIteration*stridez); //one line up in y direction, one position back in x
				v2 = *(src  -stridey + zIteration*stridez); //one line up in y direction
				v3 = *(src+1-stridey + zIteration*stridez); //one line up in y direction, one position forward in x
				v4 = *(src+2-stridey + zIteration*stridez);

				v5 = *(src-1 + zIteration*stridez); //same line in y
				v6 = *(src   + zIteration*stridez); //...
				v7 = *(src+1 + zIteration*stridez);
				v8 = *(src+2 + zIteration*stridez);

				v9  = *(src-1+stridey + zIteration*stridez);
				v10 = *(src  +stridey + zIteration*stridez);
				v11 = *(src+1+stridey + zIteration*stridez);
				v12 = *(src+2+stridey + zIteration*stridez);

				v13 = *(src-1+2*stridey + zIteration*stridez);
				v14= *(src  +2*stridey + zIteration*stridez);
				v15 = *(src+1+2*stridey + zIteration*stridez);
				v16 = *(src+2+2*stridey + zIteration*stridez);


                //printf("L-1: %.1f %.1f %.1f %.1f\\nL 0: %.1f %.1f %.1f %.1f\\nL 1: %.1f %.1f %.1f %.1f\\nL 2: %.1f %.1f %.1f %.1f\\n", v1, v2,v3,v4,v5,v6,v7,v8,v9,v10,v11,v12, v13, v14, v15, v16);


				/*calculate spline value for line 1 2 3 4 above pixel of interest */
                D1 = CSPL_INT(f1,f2,f3,f4,v1,v2,v3,v4);
                D2 = CSPL_INT(f4,f3,f2,f1,v1,v2,v3,v4);
                line1 = CSPL_CALC(v2,v3,D1,D2,xoffseth);

				/*calculate spline value for line 5 6 7 8 above pixel of interest */
                D1 = CSPL_INT(f1,f2,f3,f4,v5,v6,v7,v8);
                D2 = CSPL_INT(f4,f3,f2,f1,v5,v6,v7,v8);
                line2 = CSPL_CALC(v6,v7,D1,D2,xoffseth);

				/*calculate spline value for line 9 10 11 12 above pixel of interest */
                D1 = CSPL_INT(f1,f2,f3,f4,v9,v10,v11,v12);
                D2 = CSPL_INT(f4,f3,f2,f1,v9,v10,v11,v12);
                line3 = CSPL_CALC(v10,v11,D1,D2,xoffseth);

				/*calculate spline value for line 13 14 15 16 above pixel of interest */
                D1 = CSPL_INT(f1,f2,f3,f4,v13,v14,v15,v16);
                D2 = CSPL_INT(f4,f3,f2,f1,v13,v14,v15,v16);
                line4 = CSPL_CALC(v14,v15,D1,D2,xoffseth);





				/*finaly, calculate spline into y direction and save into value[z]*/
                D1 = CSPL_INT(f1,f2,f3,f4,line1,line2,line3,line4);
                D2 = CSPL_INT(f4,f3,f2,f1,line1,line2,line3,line4);
                layerValues[zIteration + layerOffset] = CSPL_CALC(line2,line3,D1,D2,yoffseth);

                //printf("GPU %d - %.1f %.1f %.1f - %.1f %.1f %.1f %.1f - %.1f %.1f %.1f %.1f - %.1f\\n", zIteration, x, y, z, v1, v2, v3, v4, line1, line2, line3, line4, layerValues[zIteration + layerOffset]);

			}

			if (is3D) {
				/*calculate spline value for z direction*/
				D1 = CSPL_INT(f1,f2,f3,f4,layerValues[0],layerValues[1],layerValues[2],layerValues[3]);
				D2 = CSPL_INT(f4,f3,f2,f1,layerValues[0],layerValues[1],layerValues[2],layerValues[3]);

				D1 = CSPL_CALC(layerValues[1],layerValues[2],D1,D2,zoffseth);

				return D1;

			}else {
				return layerValues[0];
			}
		}
	}
    return defaultval;
}


__device__ bool are_same(float a, float b) {

    float epsilon = 0.00001;
    return fabs(float(a) - float(b)) < epsilon;
}

__device__ bool is_identity(float *P){
    if (are_same(P[ 0], 1.) && are_same(P[ 1], 0.) && are_same(P[ 2], 0.) && are_same(P[ 3], 0.) &&
        are_same(P[ 4], 0.) && are_same(P[ 5], 1.) && are_same(P[ 6], 0.) && are_same(P[ 7], 0.) &&
        are_same(P[ 8], 0.) && are_same(P[ 9], 0.) && are_same(P[10], 1.) && are_same(P[11], 0.)) {
            return (1>0);
    } else {
        return (0 > 1);
    }
}



extern "C" __global__ 
void transformSpline(float *src, float *dst, float *P, bool is_affinity, float defaultval, float valinf, int * dim_src, int * dim_dst) {


    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;                                                                                                                                        
    int i = (blockIdx.x*(blockSize) + tid);
    
    if (i < (dim_dst[0]*dim_dst[1]*dim_dst[2])){
        int x,y,z;

        z = i  / (dim_dst[0] * dim_dst[1]);
        y = (i % (dim_dst[0] * dim_dst[1])) / dim_dst[0] ; 
        x = (i % (dim_dst[0] * dim_dst[1])) % dim_dst[0] ;


        if (is_identity(P)) {
            dst[i] = src[x + y * dim_src[0] + z * dim_src[0] * dim_src[1]];
        } else {

            float dx, dy, dz;

            dx = P[0]*x+P[1]*y+P[2]*z+P[3];
            dy = P[4]*x+P[5]*y+P[6]*z+P[7];
            dz = P[8]*x+P[9]*y+P[10]*z+P[11];

            //if (x==8 && y ==8 ){printf("SIZE GPU %d %d %d %d - %.1f %.1f %.1f; \\n", i, x,y,z, dx,dy,dz);}         

            dz = interpolCubicSpline(src, dx, dy, dz, defaultval, dim_src[0], dim_src[1], dim_src[2] );

            i = z  + y * dim_dst[2] + x * dim_dst[1] * dim_dst[2];


            //i = i;

            //printf("%d %.3f\\n", i, dz);
            dst[i] = dz;
        }    
    }
}
'''