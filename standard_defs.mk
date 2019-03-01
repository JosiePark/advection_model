
CFLAGS = -C -g -traceback -O2 -heap-arrays -check all -fp-stack-check # -g -traceback -check all -fp-stack-check -mcmodel=large
#NETCDF_VERSION=4.3.3
CC = ifort

#LIBS = -L${MKL_HOME}/interfaces/fftw3xf -L/apps/netcdf/$(NETCDF_VERSION)/lib -L${MKL_HOME}/interfaces
#LINKS = -lfftw3xf_intel -lnetcdff -lnetcdf -llapack  -lblas
LINKS = -lfftw3 -lnetcdf -lnetcdff -llapack -lblas
LIB_DIR = /usr/lib/x86_64-linux-gnu
#LIBS = -L$(LIB_DIR)

#HOME_DIR = ${HOME}/WORK/CODE
HOME_DIR = /home/clustor2/ma/j/jp1115/CODE
NETCDF_DIR = $(HOME_DIR)/NETCDF
LAGR_DIR = $(HOME_DIR)/LAGR
VAR_DIR = $(HOME_DIR)/VAR
SRC_DIR = $(HOME_DIR)/TRANSPORT
CAB_DIR = $(HOME_DIR)/CABARET
SOLVER_DIR = $(HOME_DIR)/SOLVER
INPUT_DIR = $(HOME_DIR)/INPUT
BUILD_DIR = $(HOME_DIR)/BUILD
MISC_DIR = $(HOME_DIR)/MISC
PVMAPPED_DIR = $(HOME_DIR)/PV_MAPPED
DIFF_DIR = $(HOME_DIR)/DIFFUSION_MODEL
EOF_DIR = $(HOME_DIR)/EOF_MODEL
STOCH_DIR = $(HOME_DIR)/STOCHASTIC_PARAMETERS
KINEMATIC_DIR = $(HOME_DIR)/KINEMATIC_MODEL

INCLUDE_DIR = /usr/include

INCLUDES = -I$(INCLUDE_DIR) -I$(BUILD_DIR) -I$(SRC_DIR) 
#INCLUDES = -I${MKL_HOME}/include/fftw -I${MKL_HOME}/include -I/apps/netcdf/$(NETCDF_VERSION)/include -I$(BUILD_DIR) -I$(SRC_DIR) 


