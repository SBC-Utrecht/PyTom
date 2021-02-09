extern "C" {

    __device__ float linearTex3D(cudaTextureObject_t tex, float3 coord)
    {
        return tex3D<float>(tex, coord.x, coord.y, coord.z);
    }

    __device__ float cubicTex3D(cudaTextureObject_t tex, float3 coord)
    {
        // shift the coordinate from [0,extent] to [-0.5, extent-0.5]
        const float3 coord_grid = coord - 0.5f;
        const float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        float3 w0, w1, w2, w3;
        bspline_weights(fraction, w0, w1, w2, w3);

        const float3 g0 = w0 + w1;
        const float3 g1 = w2 + w3;
        const float3 h0 = (w1 / g0) - 0.5f + index;  //h0 = w1/g0 - 1, move from [-0.5, extent-0.5] to [0, extent]
        const float3 h1 = (w3 / g1) + 1.5f + index;  //h1 = w3/g1 + 1, move from [-0.5, extent-0.5] to [0, extent]

        // fetch the eight linear interpolations
        // weighting and fetching is interleaved for performance and stability reasons
        float tex000 = tex3D<float>(tex, h0.x, h0.y, h0.z);
        float tex100 = tex3D<float>(tex, h1.x, h0.y, h0.z);
        tex000 = g0.x * tex000 + g1.x * tex100;  //weigh along the x-direction
        float tex010 = tex3D<float>(tex, h0.x, h1.y, h0.z);
        float tex110 = tex3D<float>(tex, h1.x, h1.y, h0.z);
        tex010 = g0.x * tex010 + g1.x * tex110;  //weigh along the x-direction
        tex000 = g0.y * tex000 + g1.y * tex010;  //weigh along the y-direction
        float tex001 = tex3D<float>(tex, h0.x, h0.y, h1.z);
        float tex101 = tex3D<float>(tex, h1.x, h0.y, h1.z);
        tex001 = g0.x * tex001 + g1.x * tex101;  //weigh along the x-direction
        float tex011 = tex3D<float>(tex, h0.x, h1.y, h1.z);
        float tex111 = tex3D<float>(tex, h1.x, h1.y, h1.z);
        tex011 = g0.x * tex011 + g1.x * tex111;  //weigh along the x-direction
        tex001 = g0.y * tex001 + g1.y * tex011;  //weigh along the y-direction

        return (g0.z * tex000 + g1.z * tex001);  //weigh along the z-direction
    }

    __device__ float cubicTex3DSimple(cudaTextureObject_t tex, float3 coord)
    {
        // transform the coordinate from [0,extent] to [-0.5, extent-0.5]
        const float3 coord_grid = coord - 0.5f;
        float3 index = floor(coord_grid);
        const float3 fraction = coord_grid - index;
        index = index + 0.5f;  //move from [-0.5, extent-0.5] to [0, extent]

        float result = 0.0f;
        for (float z=-1; z < 2.5f; z++)  //range [-1, 2]
        {
            float bsplineZ = bspline(z-fraction.z);
            float w = index.z + z;
            for (float y=-1; y < 2.5f; y++)
            {
                float bsplineYZ = bspline(y-fraction.y) * bsplineZ;
                float v = index.y + y;
                for (float x=-1; x < 2.5f; x++)
                {
                    float bsplineXYZ = bspline(x-fraction.x) * bsplineYZ;
                    float u = index.x + x;
                    result += bsplineXYZ * tex3D<float>(tex, u, v, w);
                }
            }
        }
        return result;
    }

}

