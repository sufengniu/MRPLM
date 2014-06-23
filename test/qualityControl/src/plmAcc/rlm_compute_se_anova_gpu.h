#include "define_t.h"
#include "cublas_v2.h"

#include "gpu_utils.h"
#include "PLMsummarize_utils.h"

void gpu_XTX(int_t y_cols, int_t y_rows, double* wts, double* XTX);

static void alloc_host_mem(int_t y_cols, int_t y_rows);

static void free_host_mem(void);

static void alloc_device_mem(int_t y_cols, int_t y_rows);

static void free_device_mem(void);

