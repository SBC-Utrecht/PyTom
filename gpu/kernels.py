spmv_kernel_text = r"""

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
void spmvh(const unsigned int Reps, const unsigned int nRow, const unsigned int prodJd, 
                      const unsigned int sumJd, const unsigned int dim, const unsigned int *Jd,
                      const unsigned int *meshindex, const unsigned int *kindx, const float2 *udata,
                      float2 *k, float2 *res, const float2 *input)
{

    int blockSize = blockDim.x; 
    unsigned int tid = threadIdx.x;
                                                                                                                                          
    unsigned int myRow0 = blockIdx.x*(blockSize) + tid;
    unsigned int myRow = myRow0/(float)Reps;
    unsigned int nc = myRow0 - myRow*Reps;
    float2 zero;
    zero.x = 0.0;
    zero.y = 0.0;
    
    if (myRow < nRow){ 
        for (unsigned int j = 0;  j  <  prodJd; j ++){
            float2 u = zero;

            // now doing the first dimension
            unsigned int index_shift = myRow * sumJd;
            
            unsigned int J = Jd[0];
            unsigned int index =    index_shift +  meshindex[dim*j + 0];
            unsigned int col = kindx[index] ;
            float2 spdata = udata[index];
            index_shift += J; 
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

            float2 ydata=input[myRow*Reps + nc]; // kout[col];
            u.x =  spdata.x*ydata.x + spdata.y*ydata.y;
            u.y =  - spdata.y*ydata.x + spdata.x*ydata.y;
            atomic_add_float2((k + col*Reps + nc), u, (res + col*Reps + nc));

        }; // Iterate for (unsigned int j = 0;  j  <  prodJd; j ++)
    };  // if (m < nRow)
};    // End of pELL_spmvh_mCoil          
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
void sumKernel(float2 *g_idata, float *weight, float * g_odata, int n) {

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
void sumKernel(float2 *g_idata, float * g_odata, int n) {

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
            
            float ix = tr[0] * float(x) + tr[2] * float(z) + float(proj_center[0]);
            float iy = tr[1] * float(y) + float(proj_center[1]);
            
            int px = int(ix); 
            int py = int(iy);
    
            float xoff = ix - float(px); 
            float yoff = iy - float(py);
            
            if (px >= 1 && py >= 1 && px < proj_dims[0] && py < proj_dims[1]){
                float v1 = projection[(px-1)*proj_dims[1]+py-1] + (projection[px*proj_dims[1]+py-1]-projection[(px-1)*proj_dims[1]+py-1]) * xoff;
                float v2 = projection[(px-1)*proj_dims[1]+py  ] + (projection[px*proj_dims[1]+py  ]-projection[(px-1)*proj_dims[1]+py  ]) * xoff;
                float v3 = v1 + (v2-v1)*yoff;
                
                atomicAdd( reconstruction + i, v3);
            };
        };
        //i++;
    //};
};
'''