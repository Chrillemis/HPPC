#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include <chrono>
#include <cstring>
#include <cuda_runtime.h>



#define BLOCK_SIZE_X 16
#define BLOCK_SIZE_Y 16

// Precision control
#define PREC 8  // 4 for float, 8 for double

#if PREC == 4
    using real_t = float;
#elif PREC == 8
    using real_t = double;
#else
    using real_t = float;
#endif



// ============================================================================
// CUDA kernels (implement the same logic as the original step routine)
// ============================================================================

/* Example
// Step 1a: ∂η/∂t = -∂(Hu)/∂x - ∂(Hv)/∂y (continuity equation)
__global__ void kernel_update_eta(real_t* eta, const real_t* u, const real_t* v,
                                   const real_t* H_at_u, const real_t* H_at_v,
                                   const char* land_mask,
                                   real_t dx, real_t dy, real_t dt, int nx, int ny) {
    ... write the code for updating eta incl land mask but excl boundaries (why ?) ...
}*/

/**
 * Advance solution by one time step using leapfrog in time:
 * 1. Update η using current u, v (continuity equation)
 * 2. Update u, v using new η (momentum equations)
 */

 /* Need to pass CUDA GPU arrays as pointers, no Array2D class available
 void step(real_t* d_eta, real_t* d_u, real_t* d_v,
          real_t g, const real_t* d_H_at_u, const real_t* d_H_at_v,
          const char* d_land_mask, const char* d_wall_u, const char* d_wall_v,
          real_t dx, real_t dy, real_t dt, int nx, int ny) {

    // declare the size of the grids and blocks that will be launched on the GPU for each kernel
    dim3 block(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    
    // we add BLOCK_SIZE_X/Y - 1 to ensure we cover the entire
    // domain even if it is not a perfect multiple of the block size
    dim3 grid((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X,
              (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y);
    // Do we need different block/grid sizes for different kernels?

    // Step 1: ∂η/∂t = -∂(Hu)/∂x - ∂(Hv)/∂y (continuity equation)
    kernel_update_eta<<<grid, block>>>(d_eta, d_u, d_v,
                                       d_H_at_u, d_H_at_v, d_land_mask,
                                       dx, dy, dt, nx, ny);

    ... rest of the kernels ...
*/

// STEP 1
__global__ void kernel_update_eta(real_t* eta, const real_t* u, const real_t* v,
    const real_t* H_at_u, const real_t* H_at_v,
    const char* land_mask,
    real_t dx, real_t dy, real_t dt, int nx, int ny) {

    int j = blockIdx.y * blockDim.y + threadIdx.y; // Row index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Column index

    if (i < nx && j < ny){
        real_t flux_x_right = H_at_u[j*(nx+1) + i+1] *u[j*(nx+1) + i+1];
        real_t flux_x_left = H_at_u[j*(nx+1) + i] *u[j*(nx+1) + i];
        real_t flux_y_top = H_at_v[(j+1)*(nx) + i] *v[(j+1)*(nx) + i];
        real_t flux_y_bottom = H_at_v[j*(nx) + i] *v[j*(nx) + i];

        real_t dflux_x = (flux_x_right - flux_x_left) / dx;
        real_t dflux_y = (flux_y_top - flux_y_bottom) / dy;

        eta[j*nx + i] = eta[j*nx + i] - dt * (dflux_x + dflux_y);
            
        // Apply land mask
        if (land_mask[j*nx + i]) {
            eta[j*nx + i] = 0.0;
        }
    }
}
        
__global__ void kernel_update_eta_boundary(real_t* eta, int nx, int ny) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // Column index

    if (idx < nx){ // Boundary check
        eta[idx] = eta[nx + idx];
        eta[nx*(ny-1) + idx] = eta[nx*(ny-2) + idx];
    }   

    if (idx < ny){
        eta[nx*idx] = eta[nx*idx + 1];
        eta[idx*nx + (nx-1)] = eta[idx*nx + (nx-2)];
    }
}

// STEP 2

__global__ void kernel_update_uv(real_t* d_eta, real_t* d_u, real_t* d_v,
    real_t g, const char* d_wall_u, const char* d_wall_v,
    real_t dx, real_t dy, real_t dt, int nx, int ny) {
    
    
    int j = blockIdx.y * blockDim.y + threadIdx.y; // Row index
    int i = blockIdx.x * blockDim.x + threadIdx.x; // Column index

    if ((0 < i) && (i < nx) && (j < ny)){ // Boundary check 
        real_t deta_dx = (d_eta[j*nx + i] - d_eta[j*(nx) + i - 1]) / dx;
        d_u[j*(nx+1) + i] = d_u[j*(nx+1)+i] - dt * g * deta_dx;
        if (d_wall_u[j*(nx+1)+i]) d_u[j*(nx+1)+i] = 0.0;

    }

    if ( (i < nx) && ( 0 < j) && (j < ny)){ // Boundary check 
        real_t deta_dy = (d_eta[j*nx + i] - d_eta[(j-1)*(nx) + i]) / dy;
        d_v[j*nx + i] = d_v[j*nx + i] - dt * g * deta_dy;
        if (d_wall_v[j*nx+i]) d_v[j*nx + i] = 0.0;

    }
}

__global__ void kernel_update_uv_boundary(real_t* d_u, real_t* d_v,
 int nx, int ny) {
    
    int idx = blockIdx.x * blockDim.x + threadIdx.x; // Column index

    if (idx < nx){ // Boundary check
        d_v[idx] = d_v[nx+idx];
        d_v[nx*(ny-1) + idx] = d_v[nx*(ny-2) + idx];
    }   

    if (idx < ny){
        d_u[(nx+1)*idx] = d_u[(nx+1)*idx + 1];
        d_u[idx*(nx+1) + (nx)] = d_u[idx*(nx+1) + (nx-1)];
    }
}


void step(real_t* d_eta, real_t* d_u, real_t* d_v,
    real_t g, const real_t* d_H_at_u, const real_t* d_H_at_v,
    const char* d_land_mask, const char* d_wall_u, const char* d_wall_v,
    real_t dx, real_t dy, real_t dt, int nx, int ny) {

    dim3 block2D(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    dim3 grid2D((nx + BLOCK_SIZE_X - 1) / BLOCK_SIZE_X,
            (ny + BLOCK_SIZE_Y - 1) / BLOCK_SIZE_Y);

    int max_dim = std::max(nx, ny);
    int block1D = 256;
    int grid1D = (max_dim + block1D - 1) / block1D;

    kernel_update_eta<<<grid2D, block2D>>>(d_eta, d_u, d_v, d_H_at_u, d_H_at_v, d_land_mask, dx, dy, dt, nx, ny);
    kernel_update_eta_boundary<<<grid1D, block1D>>>(d_eta,nx,ny);

    kernel_update_uv<<<grid2D, block2D>>>(d_eta,d_u,d_v,g,d_wall_u,d_wall_v,dx,dy,dt,nx,ny);
    kernel_update_uv_boundary<<<grid1D, block1D>>>(d_u,d_v,nx,ny);
}




void run(double duration, double output_interval,
    std::vector<real_t>& times, std::vector<real_t>& performance, 
    std::vector<Array2D<real_t>>& etas) {

    int n_steps = static_cast<int>(duration / dt);
    int output_every = output_interval > 0 ? 
                    std::max(1, static_cast<int>(output_interval / dt)) : 
                    std::max(1, n_steps / 100);

    std::cout << "Running: " << duration << "s (" << n_steps << " steps), output every " 
            << output_every * dt << "s" << std::endl;
    
    size_t size_eta  = ny * nx * sizeof(real_t);
    size_t size_u    = ny * (nx + 1) * sizeof(real_t);
    size_t size_v    = (ny + 1) * nx * sizeof(real_t);

    size_t size_mask = ny * nx * sizeof(char);
    size_t size_wall_u = ny * (nx + 1) * sizeof(char);
    size_t size_wall_v = (ny + 1) * nx * sizeof(char);

    // FIXME --- Allocate data and copy arrays to GPU ---
    real_t *d_eta, *d_u, *d_v, *d_H_at_u, *d_H_at_v;
    char *d_land_mask, *d_wall_u, *d_wall_v;

    cudaMalloc(&d_eta, size_eta);
    cudaMalloc(&d_u, size_u);
    cudaMalloc(&d_H_at_u, size_u);
    cudaMalloc(&d_v, size_v);
    cudaMalloc(&d_H_at_v, size_v);

    cudaMalloc(&d_land_mask, size_mask);
    cudaMalloc(&d_wall_u, size_wall_u);
    cudaMalloc(&d_wall_v, size_wall_v);
    


    cudaMemcpy(d_eta, eta.ptr(), size_eta, cudaMemcpyHostToDevice);
    cudaMemcpy(d_u, u.ptr(), size_u, cudaMemcpyHostToDevice);
    cudaMemcpy(d_v, v.ptr(), size_v, cudaMemcpyHostToDevice);

    cudaMemcpy(d_H_at_u, H_at_u.ptr(), size_u, cudaMemcpyHostToDevice);
    cudaMemcpy(d_H_at_v, H_at_v.ptr(), size_v, cudaMemcpyHostToDevice);

    cudaMemcpy(d_land_mask, land_mask.ptr(), size_mask, cudaMemcpyHostToDevice);
    cudaMemcpy(d_wall_u, wall_u.ptr(), size_wall_u, cudaMemcpyHostToDevice);
    cudaMemcpy(d_wall_v, wall_v.ptr(), size_wall_v, cudaMemcpyHostToDevice);

    
    // ... and more ...

    // --------------------------

    auto tlast = std::chrono::high_resolution_clock::now();
    long long nupdated = 0;

    for (int n = 0; n < n_steps; n++) {
    // FIXME Replace with call to CUDA version
    step(eta, u, v, g, H_at_u, H_at_v, land_mask, wall_u, wall_v, dx, dy, dt, nx, ny);
    nupdated += nx * ny; // increment number of updated grid points (for performance measurement)
    
    if (n % output_every == 0) {
        // We need to add a synchronise here to make sure all GPU computations are done
        cudaDeviceSynchronize();
        // FIXME copy d_eta back to CPU

        auto tnow = std::chrono::high_resolution_clock::now(); // Wall clock time since last performance measurement
        double elapsed = std::chrono::duration<double>(tnow - tlast).count();
        performance.push_back(1e9 * elapsed / nupdated); // nano-seconds per grid point update
        nupdated = 0;
        times.push_back(n * dt);
        etas.push_back(eta.copy());
        
        if (n % (10 * output_every) == 0) {
            real_t max_eta = 0;
            for (int j = 0; j < ny; j++) {
                for (int i = 0; i < nx; i++) {
                    max_eta = std::max(max_eta, std::abs(eta(j, i)));
                }
            }
            std::cout << "  Step " << n << "/" << n_steps << ", t=" << n * dt 
                        << "s, max|η|=" << max_eta << "m" << std::endl;
        }
        tlast = tnow;
    }
    }

    if (nupdated > 0) {
    // FIXME -- Same as above
    auto tnow = std::chrono::high_resolution_clock::now();
    double elapsed = std::chrono::duration<double>(tnow - tlast).count();
    performance.push_back(1e9 * elapsed / nupdated);
    times.push_back(n_steps * dt);
    etas.push_back(eta.copy());
    }

    // FIXME --- Free GPU memory ---
    cudaFree(d_eta);
    // ... and more ...
    // -----------------------

    std::cout << "Simulation complete!" << std::endl;
}

